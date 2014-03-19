function [qrs_amp_raw,qrs_i_raw] = tompkins_decision(SignalB,SignalI,Fs)

% initializations
qrs_c =[]; %amplitude of R
qrs_i =[]; %index
nois_c =[];
nois_i =[];
NOISE_LEV = 0;
skip = 0; % becomes one when a T wave is detected
not_nois = 0; % it is not noise when not_nois = 1
selected_RR =[]; % Selected RR intervals
m_selected_RR = 0;
mean_RR = 0;
qrs_i_raw =[];
qrs_amp_raw=[];
SIG_LEV = 0;
ser_back = 0;
SIG_LEV1 = 0; % Signal level in Bandpassed filter
NOISE_LEV1 = 0; % Noise level in Bandpassed filter

%% Fiducial Mark 
% Note : a minimum distance of 40 samples is considered between each R wave
% since in physiological point of view no RR wave can occur in less than
% 200 msec distance
[pks,locs] = findpeaks(SignalI,'MINPEAKDISTANCE',0.2*Fs);

%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(SignalI(1:2*Fs))*1/4; % 0.25 of the max amplitude 
THR_NOISE = mean(SignalI(1:2*Fs))*1/2; % 0.5 of the mean signal is considered to be noise

%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
THR_SIG1 = max(SignalB(1:2*Fs))*1/4; % 0.25 of the max amplitude 
THR_NOISE1 = mean(SignalB(1:2*Fs))*1/2;

%% Thresholding and online decision rule
for i = 1 : length(pks)
    
    if locs(i)-round(0.100*Fs)>= 1 && locs(i)<= length(SignalB)
        [y_i, x_i] = max(SignalB(locs(i)-30:locs(i)));
    elseif i == 1
        [y_i, x_i] = max(SignalB(1:locs(i)));
        ser_back = 1;
    elseif locs(i)>= length(SignalB)
        [y_i, x_i] = max(SignalB(locs(i)-30:end));
    end
    
    % update the heart_rate
    if length(qrs_c) >= 9
        diffRR = diff(qrs_i(end-8:end)); %calculate RR interval
        mean_RR = mean(diffRR); % calculate the mean of 8 previous R waves interval
        temp = diffRR(diffRR >= 0.92*mean_RR); %find the most regular beats
        temp = temp(temp <= 1.16*mean_RR);
        if isempty(temp)
            % lower down thresholds to detect better in MVI
            THR_SIG = 0.5*(THR_SIG);
            THR_NOISE = 0.5*(THR_NOISE);
            % lower down thresholds to detect better in Bandpass filtered 
            THR_SIG1 = 0.5*(THR_SIG1);
            THR_NOISE1 = 0.5*(THR_NOISE1);
        end
        selected_RR =[selected_RR temp];
        % keep a track of the most regular RR intervals
        if  length(selected_RR) >= 8
            m_selected_RR = mean(selected_RR);
            if length(selected_RR)-8 > 0
                selected_RR(1:length(selected_RR)-8)=[];
            end
        end
    end
    
    % find noise and QRS peaks
    if pks(i) >= THR_SIG
    
        % if a QRS candidate occurs within 360ms of the previous QRS
        % ,the algorithm determines if its T wave or QRS
        if length(qrs_c) >= 3
            if (locs(i)-qrs_i(end)) <= round(0.3600*Fs)
                Slope1 = mean(diff(SignalI(locs(i)-round(0.050*Fs):locs(i)))); %mean slope of the waveform at that position
                Slope2 = mean(diff(SignalI(qrs_i(end)-round(0.050*Fs):qrs_i(end)))); %mean slope of previous R wave
                if abs(Slope1) <= abs(0.5*(Slope2)) || (locs(i)-qrs_i(end)) <= round(0.4*test_m)  % slope less then 0.5 of previous R
                    nois_c = [nois_c pks(i)];
                    nois_i = [nois_i locs(i)];
                    skip = 1; % T wave identification
                else
                    skip = 0;
                end
            end
        end
        
        if skip == 0  % skip is 1 when a T wave is detected
            qrs_c = [qrs_c pks(i)];
            qrs_i = [qrs_i locs(i)];

            % bandpass filter check threshold
            if y_i >= THR_SIG1
                if ser_back 
                    qrs_i_raw = [qrs_i_raw x_i];  % save index of bandpass 
                else
                    qrs_i_raw = [qrs_i_raw locs(i)-30 + (x_i - 1)];% save index of bandpass 
                end
                qrs_amp_raw =[qrs_amp_raw y_i];% save amplitude of bandpass 
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;% adjust threshold for bandpass filtered sig
            end
            
            % adjust Signal level
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;
        end
        
    elseif THR_NOISE <= pks(i) && pks(i) < THR_SIG
        
        % calculate the mean of the last 8 R waves to make sure that QRS is not
        % missing(If no R detected , trigger a search back) 1.66*mean

        if m_selected_RR
            test_m = m_selected_RR; %if the regular RR availabe use it
        elseif mean_RR && m_selected_RR == 0
            test_m = mean_RR;
        else
            test_m = 0;
        end
        
        if not_nois == 0 %it is not noise when no_nois =1
            nois_c = [nois_c pks(i)];
            nois_i = [nois_i locs(i)];
            
            if THR_NOISE1 <= y_i %&&  y_i < THR_SIG1
                NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
            end
            
            %adjust Noise level
            NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
        end
    end
    
    % adjust the threshold with SNR
    if NOISE_LEV ~= 0 && SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    % adjust the threshold with SNR for bandpassed signal
    if NOISE_LEV1 ~= 0 && SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    skip = 0; %reset parameters
    not_nois = 0; %reset parameters
    ser_back = 0;  %reset bandpass param
end
