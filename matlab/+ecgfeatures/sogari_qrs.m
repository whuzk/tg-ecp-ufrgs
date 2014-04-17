function [QRSi,QRS2,THRs,THRn,RRm] = sogari_qrs(Signal,Fs)

%% initializations
N = length(Signal);             % length of signal
TTrain = 2*Fs;                  % length of training period
Tref = round(0.20*Fs);          % length of refractory period
TQtol = round(0.16*Fs);         % tolerance for comparing QRS waves
TTtol = round(0.36*Fs);         % tolerance for T wave identification
THRs = zeros(N,1);              % buffer for the Signal Threshold history
THRn = zeros(N,1);              % buffer for the Noise Threshold history
RRm = zeros(N,1);               % buffer for the RR mean history
THR_SIG = 0;                    % Signal threshold
SIG_LEV = 0;                    % Signal level
NOISE_LEV = 0;                  % Noise level
QRSi = zeros(N,1);              % buffer for the R wave indices
QRS2 = zeros(N,1);              % buffer for QRS detected in searchback
qrs_count = 0;                  % count of QRS complex in main QRS buffer
qrs_count2 = 0;                 % current position in second QRS buffer
ser_back_i = 0;                 % index of searchback starting point
last_peak_i = 0;                % index of the last peak detected
last_peak_amp = 0;              % amplitude of the last peak detected
last_qrs_i = -Inf;              % index of last detected QRS complex
qrs_amp = 0;                    % average qrs amplitude
rr_mean = Fs;                   % running average of RR intervals
rr_half = round(0.5*rr_mean);   % half of RR average
rr_miss = round(1.8*rr_mean);   % interval for qrs to be assumed as missed
signal_rising = false;          % flag to indicate a rise in the signal
lev_est_ratio = 0.125;          % ratio of signal/noise level estimation
sback_est_ratio = 0.25;         % alternative ratio for searchback
rr_est_ratio = 0.2;             % ratio of RR interval estimation
amp_est_ratio = 0.2;            % ratio of QRS amplitude estimation

%% algorithm
for i = 2:N-1
    
    % get current sample
    a = Signal(i);
    
    % initialize iteration flags
    peak_detected = false;
    qrs_detected = false;
    
    %% peak detection
    if a > last_peak_amp
        
        % signalize positive slope
        signal_rising = true;
        
        % update peak info
        last_peak_amp = a;
        last_peak_i = i;
        
    elseif ~signal_rising
        
        % update current amplitude
        last_peak_amp = a;
        
    elseif a < 0.5*last_peak_amp
        
        % signalize peak detection
        peak_detected = true;
        peak_amp = last_peak_amp;
        peak_i = last_peak_i;
        
        % reset state
        signal_rising = false;
        last_peak_amp = a;
        
    end
    
    %% qrs detection and estimation of signal/noise levels
    if peak_detected && peak_i-last_qrs_i > Tref% && peak_amp > 0.5*THR_SIG
        
        % decrease estimation ratio for times beyond the training period
        if i == TTrain
            lev_est_ratio = lev_est_ratio/2;
        end
        
        % adjust levels
        if peak_amp >= THR_SIG
            SIG_LEV = estimate(SIG_LEV,lev_est_ratio,peak_amp);
        else
            NOISE_LEV = estimate(NOISE_LEV,lev_est_ratio,peak_amp);
        end
        
        % assert qrs
    	qrs_detected = (i > TTrain) && (peak_amp >= THR_SIG) && ...
            (peak_i-last_qrs_i > max(TTtol,rr_half) || ...
            ~istwave(Signal, peak_i, last_qrs_i, TQtol));
        
    end
    
    %% search back
    if qrs_detected || ser_back_i == 0 || (signal_rising && a >= THR_SIG)
        
        % do nothing if:
        %   1. a qrs complex was detected in the current iteration
        %   2. no previous qrs were detected in previous iterations
        %   3. a new qrs might be detected in the following iteretions
        
    elseif i-ser_back_i >= 2*rr_miss
        
        SIG_LEV = 0.5*SIG_LEV;
        NOISE_LEV = 0.5*NOISE_LEV;
        ser_back_i = ser_back_i + round(rr_mean);
        %{
        if THR_SIG >= 0.01*qrs_amp
            % adjust signal and noise levels
            SIG_LEV = 0.99*SIG_LEV;
            NOISE_LEV = 0.99*NOISE_LEV;
        else
            % postpone searchback
            ser_back_i = ser_back_i + round(rr_mean);
        end
        %}
    elseif i-ser_back_i >= rr_miss
        
        % search back and locate the max in this interval
        center = ser_back_i + round(rr_mean);
        [new_amp, new_i] = findmax(Signal, center, rr_half);
        
        % check if candidate peak is from qrs
        if new_amp >= 0.5*THR_SIG
            % skip when a T wave is detected
            if new_i-last_qrs_i > max(TTtol,rr_half) || ...
                ~istwave(Signal, new_i, last_qrs_i, TQtol)
                % update peak info
                peak_i = new_i;
                peak_amp = new_amp;
                % increment qrs count
                qrs_count2 = qrs_count2 + 1;
                % save index of integrated signal
                QRS2(qrs_count2) = new_i;
                % signalize qrs detection
                qrs_detected = true;
                % adjust signal levels
                SIG_LEV = estimate(SIG_LEV,sback_est_ratio,new_amp);
            end
        end
    end
    
    
    %% update QRS info and RR-interval
    if qrs_detected
        
        % increment qrs count
        qrs_count = qrs_count + 1;
        % save index of integrated signal
        QRSi(qrs_count) = peak_i;
        
        % update RR interval
        if last_qrs_i > 0
            new_rr = peak_i - last_qrs_i;
            if 0.6*rr_mean < new_rr && new_rr < 1.4*rr_mean
                rr_mean = estimate(rr_mean,rr_est_ratio,new_rr);
                rr_mean = limit(rr_mean,0.2*Fs,2*Fs);
            end
        end
        
        % update QRS index and amplitude
        qrs_amp = estimate(qrs_amp,amp_est_ratio,peak_amp);
        last_qrs_i = peak_i;
        ser_back_i = peak_i;
        
        % calculate RR related limits
        rr_miss = round(1.8*rr_mean);
        rr_half = round(0.5*rr_mean);
    end
    
    %% update threshold
    THR_SIG = NOISE_LEV + 0.25*abs(SIG_LEV - NOISE_LEV);
    
    %% history
    THRs(i) = THR_SIG;
    THRn(i) = 0.5*THR_SIG;
    RRm(i) = rr_mean;
end

% trim the output vectors
QRSi(qrs_count+1:end) = [];
QRS2(qrs_count2+1:end) = [];


function Result = estimate(x,r,v)
Result = (1-r)*x + r*v;

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