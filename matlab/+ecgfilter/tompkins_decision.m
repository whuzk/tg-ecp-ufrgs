function [QRSb,QRSi,THRa] = tompkins_decision(SignalB,SignalI,DelayI,Fs)

%% initializations
N1 = length(SignalB);       % length of bandpassed signal
N2 = length(SignalI);       % length of MWI signal
RefPeriod = round(0.2*Fs);  % length of refractory period
TrainPeriod = round(2*Fs);  % length of training period

QRSb = zeros(N1,1);         % buffer for the R wave indices
QRSi = zeros(N1,1);         % buffer for the R wave indices
THRa = zeros(N2,1);         % buffer for the Threshold history
RR_int = NaN(1,8);          % buffer for the last RR intervals
sel_RR_int = NaN(1,8);      % buffer for the selected RR intervals

detection = false;          % flag to indicate the detection phase
qrs_updated = false;        % flag to indicate when a new qrs is detected
rr_mean = 0;                % running average of latest RR intervals
rr_count = 0;               % count of selected RR intervals in RR buffer
sel_rr_count = 0;           % count of selected RR intervals in RR buffer
qrs_count = 0;              % count of QRS complex in output QRS buffer
last_peak_i = 0;            % index of last peak
last_qrs_i = 0;             % index of last QRS complex
new_rr = 0;                 % new RR interval after QRS detection
next_limit = 0;             % limit before which peaks must be discarded
                            % due to the refractory period

SIG_LEV1 = 0;               % Signal level in MWI
NOISE_LEV1 = 0;             % Noise level in MWI
SIG_LEV2 = 0;               % Signal level in Bandpass
NOISE_LEV2 = 0;             % Noise level in Bandpass

THR_SIG1 = 0;               % Signal threshold in MWI
THR_NOISE1 = 0;             % Noise threshold in MWI
THR_SIG2 = 0;               % Signal threshold in Bandpass
THR_NOISE2 = 0;             % Noise threshold in Bandpass

%% Find peaks, adapt thresholds and perform online decision rule
[pks,locs] = findpeaks(SignalI);

for i = 1:length(pks)
    
    % get current peak location
    peak_a = pks(i);
    peak_i = locs(i);
    
    % check if we are in the refractory period
    if peak_i <= next_limit
        continue;
    end
    
    % update threshold history
    THRa(last_peak_i+1:peak_i) = THR_SIG1;
    last_peak_i = peak_i;
    
    % find peak location in the bandpassed signal
    begin_i = max(0,peak_i-2*DelayI);
    end_i = min(peak_i,length(SignalB));
    [y_i, x_i] = max(SignalB(begin_i:end_i));
    
    %% check whether a new QRS was detected and update heart rates
    if qrs_updated && peak_i > TrainPeriod && new_rr > 0
        if ~detection
            % initiate detection phase
            detection = true;
        end
        
        % update RR intervals
        rr_count = mod(rr_count,length(RR_int))+1;
        RR_int(rr_count) = new_rr;
        
        % calculate running average of RR intervals
        rr_mean = nanmean(RR_int);
        is_regular = (0.92*rr_mean <= new_rr && new_rr <= 1.16*rr_mean);

        % check if new beat is regular
        if is_regular
            % update buffer of selected RR intervals
            sel_rr_count = mod(sel_rr_count,length(sel_RR_int))+1;
            sel_RR_int(sel_rr_count) = new_rr;
            
            % substitute the average
            rr_mean = nanmean(sel_RR_int);
        else
            % reduce thresholds to detect better in MVI
            THR_SIG1 = 0.5*THR_SIG1;
            THR_NOISE1 = 0.5*THR_NOISE1;
            
            % reduce thresholds to detect better in Bandpass filtered 
            THR_SIG2 = 0.5*THR_SIG2;
            THR_NOISE2 = 0.5*THR_NOISE2;
        end
    end
    
    % reset QRS detection flag
    qrs_updated = false;
    
    %% find noise and QRS peaks
    if peak_a >= THR_SIG1
    
        % reset T wave indication flag
        is_twave = false;
        
        % if a QRS candidate occurs within 360ms of the previous QRS
        % the algorithm determines if its T wave or QRS
        if detection && (peak_i-last_qrs_i) <= round(0.36*Fs)
            %mean slope of the waveform at that position
            Slope1 = mean(diff(SignalI(peak_i-round(0.05*Fs):peak_i)));
            %mean slope of previous R wave
            Slope2 = mean(diff(SignalI(last_qrs_i-round(0.05*Fs):last_qrs_i)));
            % slope less then 0.5 of previous R
            is_twave = abs(Slope1) <= abs(0.5*(Slope2)) || ...
                (peak_i-last_qrs_i) <= round(0.4*rr_mean);
        end
        
        % skip when a T wave is detected
        if ~is_twave
            % update last QRS peak and RR interval
            if last_qrs_i > 0
                new_rr = peak_i - last_qrs_i;
            end
            last_qrs_i = peak_i;
            qrs_updated = true;
            
            % bandpass filter check threshold
            if y_i >= THR_SIG2
                if detection
                    % increment qrs count
                    qrs_count = qrs_count + 1;
                    % save index of bandpass
                    QRSb(qrs_count) = begin_i + x_i - 1;
                    % save index of MWI 
                    QRSi(qrs_count) = peak_i;
                    % set next limit
                    next_limit = peak_i + RefPeriod;
                end
                % adjust threshold for bandpass signal
                SIG_LEV2 = 0.125*y_i + 0.875*SIG_LEV2;
            end
            
            % adjust Signal level
            SIG_LEV1 = 0.125*peak_a + 0.875*SIG_LEV1;
        end
        
    elseif THR_NOISE1 <= peak_a && peak_a < THR_SIG1
        
        % calculate the mean of the last 8 R waves to make sure that QRS is
        % not missing (if no R detected, trigger a search back) 1.66*mean
        if detection && (peak_i-last_qrs_i) >= round(1.66*rr_mean)
            % search back and locate the max in this interval
            peak_a = max(SignalI(last_qrs_i + RefPeriod:peak_i));
            if peak_a >= THR_NOISE1
                % update last QRS peak and RR interval
                if last_qrs_i > 0
                    new_rr = peak_i - last_qrs_i;
                end
                last_qrs_i = peak_i;
                qrs_updated = true;

                % take care of bandpass signal threshold
                if y_i >= THR_NOISE2
                    % increment qrs count
                    qrs_count = qrs_count + 1;
                    % save index of bandpass 
                    QRSb(qrs_count) = begin_i + x_i - 1;
                    % save index of MWI 
                    QRSi(qrs_count) = peak_i;
                    % set next limit
                    next_limit = peak_i + RefPeriod;
                    % when found with the second threshold
                    SIG_LEV2 = 0.25*y_i + 0.75*SIG_LEV2;
                end

                %when found with the second threshold
                SIG_LEV1 = 0.25*peak_a + 0.75*SIG_LEV1;
            end
        else
            %adjust noise level 1
            NOISE_LEV1 = 0.125*peak_a + 0.875*NOISE_LEV1;
            %adjust noise level 2
            if THR_NOISE2 <= y_i %&&  y_i < THR_SIG1
                NOISE_LEV2 = 0.125*y_i + 0.875*NOISE_LEV2;
            end
        end
    end
    
    %% adjust thresholds
    % adjust the threshold with SNR for MWI signal
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_NOISE1 = 0.5*THR_SIG1;
    
    % adjust the threshold with SNR for bandpassed signal
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    THR_NOISE2 = 0.5*THR_SIG2;
end

% complete the threshold history
THRa(last_peak_i+1:end) = THR_SIG1;

% trim the output vectors
QRSb(qrs_count+1:end) = [];
QRSi(qrs_count+1:end) = [];