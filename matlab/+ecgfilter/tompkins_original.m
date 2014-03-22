function [QRSi,QRS2,THRs,THRn,THRs2,THRn2] = tompkins_original(SignalB,SignalI,DelayI,Fs)

%% initializations
N = length(SignalI);    % length of MWI signal
N2 = length(SignalB);   % length of bandpassed signal
Tref = round(0.2*Fs);   % length of refractory period
Ttrain = round(2*Fs);   % length of training period
Tqrs = round(0.05*Fs);  % half the length of a QRS wave
TTtol = round(0.36*Fs); % tolerance for T wave identification

QRSi = zeros(N,1);      % buffer for the R wave indices
QRS2 = zeros(N,1);      % buffer for the R wave indices
THRs = zeros(N,1);      % buffer for the Signal Threshold history
THRn = zeros(N,1);      % buffer for the Noise Threshold history
THRs2 = zeros(N2,1);    % buffer for the Signal Threshold history
THRn2 = zeros(N2,1);    % buffer for the Noise Threshold history
RR_int = NaN(1,8);      % buffer for the last RR intervals
sel_RR_int = NaN(1,8);  % buffer for the selected RR intervals

acthung = false;        % flag to indicate the active search
detection = false;      % flag to indicate the detection phase
qrs_updated = false;    % flag to indicate when a new qrs is detected
rr_mean = 0;            % running average of latest RR intervals
rr_count = 0;           % count of selected RR intervals in RR buffer
sel_rr_count = 0;       % count of selected RR intervals in RR buffer
qrs_count = 0;          % count of QRS complex in output QRS buffer
qrs_count2 = 0;         % count of QRS complex in output QRS buffer
last_peak_i = 0;        % index of last identified peak
cur_i = 0;              % index of current candidate peak
cur_a = 0;              % amplitude of current candidate peak in MWI
cur_y = 0;              % amplitude of current candidate peak in BP
last_qrs_i = 0;         % index of last QRS complex
new_rr = 0;             % new RR interval after QRS detection
ref_count = 0;          % counter for the refractory period

SIG_LEV1 = 0;           % Signal level in MWI
NOISE_LEV1 = 0;         % Noise level in MWI
SIG_LEV2 = 0;           % Signal level in Bandpass
NOISE_LEV2 = 0;         % Noise level in Bandpass

THR_SIG1 = 0;           % Signal threshold in MWI
THR_NOISE1 = 0;         % Noise threshold in MWI
THR_SIG2 = 0;           % Signal threshold in Bandpass
THR_NOISE2 = 0;         % Noise threshold in Bandpass

%% algorithm
for i = 2:length(SignalI)-1
    
    % get current and neighbour samples
    a = SignalI(i);
    b = SignalI(i-1);
    c = SignalI(i+1);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif acthung && a < 0.5*cur_a
    
        if detection
            % reset T wave indication flag
            is_twave = false;
            
            % if a QRS candidate occurs within 360ms of the previous QRS
            % the algorithm determines if it is T wave or QRS
            if cur_i-last_qrs_i < TTtol
                %mean slope of the waveform at that position
                range1 = cur_i-Tqrs:cur_i+Tqrs;
                slope1 = max(diff(SignalI(range1)));
                %mean slope of previous R wave
                range2 = last_qrs_i-Tqrs:last_qrs_i+Tqrs;
                slope2 = max(diff(SignalI(range2)));
                % slope less then 0.5 of previous R
                is_twave = abs(slope1) < abs(0.5*slope2);
            end
            
            % skip when a T wave is detected
            if ~is_twave && cur_y >= THR_SIG2
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of MWI
                QRSi(qrs_count) = cur_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-cur_i));
            end
        end
        
        % update last QRS peak and RR interval
        if last_qrs_i > 0
            new_rr = cur_i - last_qrs_i;
        end
        last_qrs_i = cur_i;
        qrs_updated = true;
        acthung = false;
        
    elseif b-a < 0 && a-c >= 0
    
        % update threshold history
        THRs(last_peak_i+1:i) = THR_SIG1;
        THRn(last_peak_i+1:i) = THR_NOISE1;
        THRs2(last_peak_i+1:i) = THR_SIG2;
        THRn2(last_peak_i+1:i) = THR_NOISE2;
        last_peak_i = i;

        % find peak in the bandpassed signal
        begin_i = max(1,i-2*DelayI);
        end_i = min(i,length(SignalB));
        y = max(SignalB(begin_i:end_i));
        
        %% find noise and QRS peaks
        if a >= THR_SIG1
            
            % check if peak is higher
            if ~acthung || a > cur_a
                cur_a = a;
                cur_i = i;
                cur_y = y;
                acthung = true;
            end
            
            % bandpass filter check threshold
            if y >= THR_SIG2
                % adjust threshold for bandpass signal
                SIG_LEV2 = 0.125*y + 0.875*SIG_LEV2;
            end
            % adjust Signal level
            SIG_LEV1 = 0.125*a + 0.875*SIG_LEV1;
            
        elseif a >= THR_NOISE1
            
            %adjust noise level 1
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
            %adjust noise level 2
            if THR_NOISE2 <= y && y < THR_SIG2
                NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
            end
            
        end
    
    elseif detection && (a < THR_NOISE1) && (i-last_qrs_i) >= round(1.66*rr_mean)
        
        % search back and locate the max in this interval
        begin_i = last_qrs_i + Tref;
        [new_a, x] = max(SignalI(begin_i:i));
        new_i = begin_i + x - 1;
        
        if new_a >= THR_NOISE1 %&& new_a < THR_SIG1
            
            % update last QRS peak and RR interval
            if last_qrs_i > 0
                new_rr = new_i - last_qrs_i;
            end
            last_qrs_i = new_i;
            qrs_updated = true;
            
            % find peak in the bandpassed signal
            begin_i = max(1,new_i-2*DelayI);
            end_i = min(new_i,length(SignalB));
            new_y = max(SignalB(begin_i:end_i));
            
            % take care of bandpass signal threshold
            if new_y >= THR_NOISE2
                % increment qrs count
                qrs_count2 = qrs_count2 + 1;
                % save index of MWI 
                QRS2(qrs_count2) = new_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-new_i));
                % when found with the second threshold
                SIG_LEV2 = 0.25*new_y + 0.75*SIG_LEV2;
            end
            
            %when found with the second threshold
            SIG_LEV1 = 0.25*new_a + 0.75*SIG_LEV1;
        end
    end
    
    %% check whether a new QRS was detected and update heart rates
    if qrs_updated && i > Ttrain && new_rr > 0
        if ~detection
            % initiate detection phase
            detection = true;
        end
        
        % update RR intervals calculate running average
        rr_count = mod(rr_count,length(RR_int))+1;
        RR_int(rr_count) = new_rr;
        rr_mean = nanmean(RR_int);
        
        % check if new beat is regular
        if new_rr >= 0.92*rr_mean && new_rr <= 1.16*rr_mean
            % update buffer of selected RR intervals
            sel_rr_count = mod(sel_rr_count,length(sel_RR_int))+1;
            sel_RR_int(sel_rr_count) = new_rr;
        end
        
        % check if heart rate is regular
        irregular = RR_int < 0.92*rr_mean | RR_int > 1.16*rr_mean;
        if ~isempty(find(irregular, 1))
            % substitute the average
            rr_mean = nanmean(sel_RR_int);
            
            % reduce thresholds to detect better in MVI
            SIG_LEV1 = 0.5*SIG_LEV1;
            NOISE_LEV1 = 0.5*NOISE_LEV1;
            
            % reduce thresholds to detect better in Bandpass filtered 
            SIG_LEV2 = 0.5*SIG_LEV2;
            NOISE_LEV2 = 0.5*NOISE_LEV2;
        end
    end
    
    % reset QRS detection flag
    qrs_updated = false;
    
    %% adjust thresholds
    % adjust the threshold with SNR for MWI signal
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_NOISE1 = 0.5*THR_SIG1;

    % adjust the threshold with SNR for bandpassed signal
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    THR_NOISE2 = 0.5*THR_SIG2;
end

% complete the threshold history
THRs(last_peak_i+1:end) = THR_SIG1;
THRn(last_peak_i+1:end) = THR_NOISE1;
THRs2(last_peak_i+1:end) = THR_SIG2;
THRn2(last_peak_i+1:end) = THR_NOISE2;

% trim the output vectors
QRSi(qrs_count+1:end) = [];
QRS2(qrs_count2+1:end) = [];