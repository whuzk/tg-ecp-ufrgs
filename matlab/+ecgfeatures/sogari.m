function [QRSi,QRS2,THRs,THRn,RRm] = sogari(Signal,Fs)

%% global initializations
D = [0; diff(Signal)];      % derivative of the integrated signal
N = length(Signal);         % length of integrated signal
Ttrain = 2*Fs;              % length of training period
Tref = round(0.20*Fs);      % length of refractory period
TQtol = round(0.16*Fs);     % tolerance for comparing QRS waves
TTtol = round(0.36*Fs);     % tolerance for T wave identification
TPtol = round(0.05*Fs);     % tolerance for peaks in filtered signal

%% initializations for phase 1
THRs = zeros(N,1);      % buffer for the Signal Threshold history
THRn = zeros(N,1);      % buffer for the Noise Threshold history
RRm = zeros(N,1);       % buffer for the RR mean history
THR_SIG = 0;            % Signal threshold in integrated signal
SIG_LEV = 0;            % Signal level in integrated signal
NOISE_LEV = 0;          % Noise level in integrated signal

%% algorithm - phase 1
for i = TPtol+1:Ttrain
    
    % check if current sample is a peak
    if D(i) > 0 && D(i+1) <= 0
        
        % get current sample
        a = Signal(i);
        
        % check if peak is from qrs
        if a >= THR_SIG
            % adjust signal levels
            SIG_LEV = 0.25*a + 0.75*SIG_LEV;
        elseif a >= 0.5*THR_SIG
            % adjust noise levels
            NOISE_LEV = 0.25*a + 0.75*NOISE_LEV;
        end
    end
    
    % update thresholds
    THR_SIG = NOISE_LEV + 0.25*abs(SIG_LEV - NOISE_LEV);
    
    % update threshold history
    THRs(i) = THR_SIG;
    THRn(i) = 0.5*THR_SIG;
end

%% initializations for phase 2
QRSi = zeros(N,1);      % buffer for the R wave indices
act_search = false;     % flag to indicate the active search
ref_count = 0;          % counter for the refractory period
qrs_count = 0;          % count of QRS complex in main QRS buffer
last_qrs_i = 0;         % index of last detected QRS complex
cur_i = 0;              % index of current candidate peak
cur_a = 0;              % amplitude of current candidate peak in integrated signal
new_rr = 0;             % new RR interval after QRS detection
rr_mean = Fs;           % running average of RR intervals
qrs_amp = 0;            % average qrs amplitude

%% algorithm - phase 2
i = Ttrain+1;
%{
while qrs_count < 3 && i < N-TPtol+1
    
    % get current sample
    a = Signal(i);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif act_search && a < 0.5*cur_a
        
        % increment qrs count
        qrs_count = qrs_count + 1;
        % save index of integrated signal
        QRSi(qrs_count) = cur_i;
        % activate refractory period
        ref_count = max(0,Tref-(i-cur_i));
        % update last QRS index and RR interval
        if last_qrs_i > 0
            new_rr = cur_i - last_qrs_i;
            rr_mean = limit(max(rr_mean,new_rr),0.2*Fs,2*Fs);
            qrs_amp = max(qrs_amp,cur_a);
        end
        last_qrs_i = cur_i;
        % adjust signal levels
        SIG_LEV = 0.125*cur_a + 0.875*SIG_LEV;
        % reset active search flag
        act_search = false;
    
    elseif D(i) > 0 && D(i+1) <= 0
        
        % check if peak is candidate to be from qrs
        if a >= THR_SIG
            % trigger or update the active search
            if ~act_search || a > cur_a
                cur_a = a;
                cur_i = i;
                act_search = true;
            end
        elseif a >= 0.5*THR_SIG
            % adjust noise levels
            NOISE_LEV = 0.125*a + 0.875*NOISE_LEV;
        end
    end
    
    % update thresholds
    THR_SIG = NOISE_LEV + 0.25*abs(SIG_LEV - NOISE_LEV);
    
    % update threshold history
    THRs(i) = THR_SIG;
    THRn(i) = 0.5*THR_SIG;
    RRm(i) = rr_mean;
    
    i = i + 1;
end
%}
%% initializations for phase 3
QRS2 = zeros(N,1);              % buffer for QRS detected in searchback
qrs_count2 = 0;                 % current position in second QRS buffer
rr_half = round(0.5*rr_mean);   % half of RR average
rr_miss = round(1.8*rr_mean);   % interval for qrs to be assumed missed
qrs_updated = false;            % flag to indicate detection of QRS
ser_back_i = i;                 % index of searchback starting point

%% algorithm - phase 3
for i = i:N-TPtol
    
    % get current sample
    a = Signal(i);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif act_search && a < 0.5*cur_a
    
        % skip when a T wave is detected
        if cur_i-last_qrs_i > max(TTtol,rr_half) || ...
            ~istwave(Signal, cur_i, last_qrs_i, TQtol)
            % increment qrs count
            qrs_count = qrs_count + 1;
            % save index of integrated signal
            QRSi(qrs_count) = cur_i;
            % activate refractory period
            ref_count = max(0,Tref-(i-cur_i));
            % update last QRS index and RR interval
            if last_qrs_i > 0
                new_rr = cur_i - last_qrs_i;
            end
            last_qrs_i = cur_i;
            qrs_updated = true;
            % adjust signal levels
            SIG_LEV = 0.125*cur_a + 0.875*SIG_LEV;
        end
        
        % reset active search flag
        act_search = false;
        
    elseif D(i) > 0 && D(i+1) <= 0
    
        % check if peak is candidate to be from qrs
        if a >= THR_SIG
            % trigger or update the active search
            if ~act_search || a > cur_a
                cur_a = a;
                cur_i = i;
                act_search = true;
                ser_back_i = i;
            end
        elseif a >= 0.5*THR_SIG
            % adjust noise levels
            NOISE_LEV = 0.125*a + 0.875*NOISE_LEV;
        end
    
    elseif ser_back_i > 0 && i-ser_back_i >= 2*rr_miss
        
        if THR_SIG >= 0.01*qrs_amp
            % adjust noise levels
            SIG_LEV = 0.99*SIG_LEV;
            NOISE_LEV = 0.99*NOISE_LEV;
        else
            % postpone searchback
            ser_back_i = ser_back_i + round(rr_mean);
        end
        
    elseif ser_back_i > 0 && i-ser_back_i >= rr_miss
        
        % search back and locate the max in this interval
        center = ser_back_i + round(rr_mean);
        [new_a, new_i] = findmax(Signal, center, rr_half);
        
        % check if candidate peak is from qrs
        if new_a >= 0.5*THR_SIG
            % skip when a T wave is detected
            if new_i-last_qrs_i > max(TTtol,rr_half) || ...
                ~istwave(Signal, new_i, last_qrs_i, TQtol)
                % update last QRS peak and RR interval
                if last_qrs_i > 0
                    new_rr = new_i - last_qrs_i;
                end
                last_qrs_i = new_i;
                qrs_updated = true;
                % increment qrs count
                qrs_count2 = qrs_count2 + 1;
                % save index of integrated signal
                QRS2(qrs_count2) = new_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-new_i));
                % adjust signal levels
                SIG_LEV = 0.25*new_a + 0.75*SIG_LEV;
            end
        end
    end
    
    % check whether a new QRS was detected
    if qrs_updated
        
        % update RR interval average
        if new_rr > 0 && 0.6*rr_mean < new_rr && new_rr < 1.4*rr_mean
            rr_mean = limit(0.8*rr_mean + 0.2*new_rr,0.2*Fs,2*Fs);
        end
        
        % update qrs amplitude
        qrs_amp = 0.8*qrs_amp + 0.2*Signal(last_qrs_i);
        
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
    THR_SIG = NOISE_LEV + 0.25*abs(SIG_LEV - NOISE_LEV);
    
    % update threshold history
    THRs(i) = THR_SIG;
    THRn(i) = 0.5*THR_SIG;
    RRm(i) = rr_mean;
end

% trim the output vectors
QRSi(qrs_count+1:end) = [];
QRS2(qrs_count2+1:end) = [];


function Result = limit(x,a,b)
Result = min(max(x,a),b);

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
% max slope less then that of previous QRS
Result = slope1 < slope2;