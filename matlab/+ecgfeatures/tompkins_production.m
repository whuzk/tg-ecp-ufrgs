function QRS = tompkins_production(SignalF,SignalI,Fs)

%% global initializations
SignalF = abs(SignalF);     % absolute value of filtered signal
D = [0; diff(SignalI)];     % derivative of the integrated signal
N = length(SignalI);        % length of integrated signal
Ttrain = 2*Fs;              % length of training period
Tref = round(0.20*Fs);      % length of refractory period
TQtol = round(0.10*Fs);     % tolerance for comparing QRS waves
TTtol = round(0.36*Fs);     % tolerance for T wave identification
TPtol = round(0.05*Fs);     % tolerance for peaks in filtered signal

%% initializations for phase 1
THR_SIG1 = 0;           % Signal threshold in integrated signal
SIG_LEV1 = 0;           % Signal level in integrated signal
NOISE_LEV1 = 0;         % Noise level in integrated signal
THR_SIG2 = 0;           % Signal threshold in filtered signal
SIG_LEV2 = 0;           % Signal level in filtered signal
NOISE_LEV2 = 0;         % Noise level in filtered signal

%% algorithm - phase 1
for i = TPtol+1:Ttrain
    
    % check if current sample is a peak
    if D(i) > 0 && D(i+1) <= 0
        
        % get current sample
        a = SignalI(i);
        
        % find peak in the filtered signal
        y = findmax(SignalF, i, TPtol);
        
        % check if peak is from qrs
        if a >= THR_SIG1 && y >= THR_SIG2
            % adjust signal levels
            SIG_LEV1 = 0.25*a + 0.75*SIG_LEV1;
            SIG_LEV2 = 0.25*y + 0.75*SIG_LEV2;
        else
            % adjust noise levels
            NOISE_LEV1 = 0.25*a + 0.75*NOISE_LEV1;
            NOISE_LEV2 = 0.25*y + 0.75*NOISE_LEV2;
        end
    end
    
    % update thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
end

%% initializations for phase 2
QRS = zeros(N,1);      % buffer for the R wave indices
act_search = false;     % flag to indicate the active search
ref_count = 0;          % counter for the refractory period
qrs_count = 0;          % count of QRS complex in main QRS buffer
last_qrs_i = 0;         % index of last detected QRS complex
cur_i = 0;              % index of current candidate peak
cur_a = 0;              % amplitude of current candidate peak in integrated signal
cur_y = 0;              % amplitude of current candidate peak in filtered signal
new_rr = 0;             % new RR interval after QRS detection

%% algorithm - phase 2
i = Ttrain+1;
while new_rr == 0 && i < N-TPtol+1
    
    % get current sample
    a = SignalI(i);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif act_search && a < 0.5*cur_a
        
        % check if candidate peak is from qrs
        if cur_y >= THR_SIG2
            % increment qrs count
            qrs_count = qrs_count + 1;
            % save index of integrated signal
            QRS(qrs_count) = cur_i;
            % activate refractory period
            ref_count = max(0,Tref-(i-cur_i));
            % update last QRS index and RR interval
            if last_qrs_i > 0
                new_rr = cur_i - last_qrs_i;
            end
            last_qrs_i = cur_i;
            % adjust signal levels
            SIG_LEV1 = 0.125*cur_a + 0.875*SIG_LEV1;
            SIG_LEV2 = 0.125*cur_y + 0.875*SIG_LEV2;
        end
        
        % reset active search flag
        act_search = false;
    
    elseif D(i) > 0 && D(i+1) <= 0
        
        % find peak in the filtered signal
        y = findmax(SignalF, i, TPtol);
        
        % check if peak is candidate to be from qrs
        if a >= THR_SIG1
            % trigger or update the active search
            if ~act_search || a > cur_a
                cur_a = a;
                cur_i = i;
                cur_y = y;
                act_search = true;
            end
        else
            % adjust noise levels
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
    end
    
    % update thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    
    i = i + 1;
end

%% initializations for phase 3
rr_mean = new_rr;               % running avereage of RR intervals
rr_half = round(0.5*rr_mean);   % half of RR average
rr_miss = round(1.8*rr_mean);  % interval for qrs to be assumed missed
qrs_updated = false;            % flag to indicate detection of QRS
ser_back_i = 0;                 % index of searchback starting point

%% algorithm - phase 3
for i = i:N-TPtol
    
    % get current sample
    a = SignalI(i);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif act_search && a < 0.5*cur_a
    
        % check if candidate peak is from qrs
        if cur_y >= THR_SIG2
            % skip when a T wave is detected
            if cur_i-last_qrs_i > TTtol || ...
                ~istwave(SignalI, cur_i, last_qrs_i, TQtol)
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of integrated signal
                QRS(qrs_count) = cur_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-cur_i));
                % update last QRS index and RR interval
                new_rr = cur_i - last_qrs_i;
                last_qrs_i = cur_i;
                qrs_updated = true;
                % adjust signal levels
                SIG_LEV1 = 0.125*cur_a + 0.875*SIG_LEV1;
                SIG_LEV2 = 0.125*cur_y + 0.875*SIG_LEV2;
            end
        end
        
        % reset active search flag
        act_search = false;
        
    elseif D(i) > 0 && D(i+1) <= 0
    
        % find peak in the filtered signal
        y = findmax(SignalF, i, TPtol);
        
        % check if peak is candidate to be from qrs
        if a >= THR_SIG1
            % trigger or update the active search
            if ~act_search || a > cur_a
                cur_a = a;
                cur_i = i;
                cur_y = y;
                act_search = true;
                ser_back_i = i;
            end
        else
            % adjust noise levels
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
        
    elseif ser_back_i > 0 && i-ser_back_i >= rr_miss
        
        % search back and locate the max in this interval
        center = ser_back_i + round(rr_mean);
        [new_a, new_i] = findmax(SignalI, center, rr_half);
        
        % find peak in the filtered signal
        new_y = findmax(SignalF, new_i, TPtol);
        
        % check if candidate peak is from qrs
        if new_a >= 0.5*THR_SIG1 && new_y >= 0.5*THR_SIG2
            % skip when a T wave is detected
            if new_i-last_qrs_i > TTtol || ...
                ~istwave(SignalI, new_i, last_qrs_i, TQtol)
                % update last QRS peak and RR interval
                new_rr = new_i - last_qrs_i;
                last_qrs_i = new_i;
                qrs_updated = true;
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of integrated signal
                QRS(qrs_count) = new_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-new_i));
            end
        end
        
        % check if qrs was found
        if qrs_updated
            % adjust signal levels
            SIG_LEV1 = 0.25*new_a + 0.75*SIG_LEV1;
            SIG_LEV2 = 0.25*new_y + 0.75*SIG_LEV2;
        elseif THR_SIG1 >= 0.05*SignalI(last_qrs_i)
            % adjust noise levels
            SIG_LEV1 = 0.975*SIG_LEV1;
            SIG_LEV2 = 0.975*SIG_LEV2;
            NOISE_LEV1 = 0.975*NOISE_LEV1;
            NOISE_LEV2 = 0.975*NOISE_LEV2;
        else
            % postpone searchback
            ser_back_i = center;
        end
        
    end
    
    % check whether a new QRS was detected
    if qrs_updated
        
        % update RR interval average
        if 0.8*rr_mean < new_rr && new_rr < 1.2*rr_mean
            rr_mean = 0.8*rr_mean + 0.2*new_rr;
        end
        
        % calculate RR miss limit
        rr_miss = round(1.8*rr_mean);
        % calculate half of RR mean
        rr_half = round(0.5*rr_mean);
        
        % enable searchback
        ser_back_i = last_qrs_i;
        % reset detection flag
        qrs_updated = false;
    end
    
    % update thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
end

% trim the output vectors
QRS(qrs_count+1:end) = [];


function [y,x] = findmax(Signal,i,l)
begin = max(1,i-l);
[y,x] = max(Signal(begin:i+l));
x = x + begin - 1;

function Result = istwave(Signal, candQrs, lastQrs, qrsLen)
% max slope of the candidate waveform
begin1 = max(1,candQrs-qrsLen);
slope1 = max(diff(Signal(begin1:candQrs)));
% max slope of previous QRS waveform
begin2 = max(1,lastQrs-qrsLen);
slope2 = max(diff(Signal(begin2:lastQrs)));
% slope less then 0.5 of previous QRS
Result = slope1 < 0.5*slope2;