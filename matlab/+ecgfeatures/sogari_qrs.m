function [QRSi,QRS2,THRs,THRn,RRm] = sogari_qrs(Signal,Fs)

%% initializations
N = length(Signal);             % length of signal
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
rr_miss = round(1.8*rr_mean);   % interval for qrs to be assumed as missed
signal_rising = false;          % flag to indicate a rise in the signal
lev_est_ratio = 0.125;          % ratio of signal/noise level estimation

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
    if peak_detected && peak_i-last_qrs_i > round(0.2*Fs)
        
        % decrease estimation ratio for times beyond the training period
        if i == 2*Fs
            lev_est_ratio = lev_est_ratio/2;
        end
        
        % adjust levels
        if peak_amp >= THR_SIG
            SIG_LEV = estimate(SIG_LEV,lev_est_ratio,peak_amp);
        else
            NOISE_LEV = estimate(NOISE_LEV,lev_est_ratio,peak_amp);
        end
        
        % assert qrs
        if (i > 2*Fs) && (peak_amp >= THR_SIG) && ~istwave(Signal, peak_i, last_qrs_i, Fs)
            qrs_detected = true;
        end
    end
    
    %% search back
    if qrs_detected || (signal_rising && a >= THR_SIG)
        
        % do nothing if:
        %   1. a qrs complex was detected in the current iteration
        %   2. a new qrs might be detected in the following iterations
        
    elseif (ser_back_i > 0) && i-ser_back_i >= rr_miss
        
        % search back and locate the max in this interval
        center = ser_back_i + round(rr_mean);
        [peak_amp, peak_i] = findmax(Signal, center, round(0.5*rr_mean));
        
        % check if candidate peak is from qrs
        if peak_amp >= 0.5*THR_SIG && ~istwave(Signal, peak_i, last_qrs_i, Fs)
            % increment qrs count
            qrs_count2 = qrs_count2 + 1;
            % save index of integrated signal
            QRS2(qrs_count2) = peak_i;
            % signalize qrs detection
            qrs_detected = true;
            % adjust signal levels
            SIG_LEV = estimate(SIG_LEV,0.25,peak_amp);
        else
            % reduce levels by half
            SIG_LEV = 0.5*SIG_LEV;
            NOISE_LEV = 0.5*NOISE_LEV;
            % postpone searchback
            ser_back_i = center;
        end
    end
    
    
    %% update QRS info and RR-interval
    if qrs_detected
        
        % push new QRS index to output buffer
        qrs_count = qrs_count + 1;
        QRSi(qrs_count) = peak_i;
        
        % update QRS amplitude estimation
        qrs_amp = estimate(qrs_amp,0.2,peak_amp);
        
        % update RR interval
        if last_qrs_i > 0
            new_rr = peak_i - last_qrs_i;
            if 0.6*rr_mean < new_rr && new_rr < 1.4*rr_mean
                rr_mean = estimate(rr_mean,0.2,new_rr);
                rr_mean = limit(rr_mean,0.2*Fs,2*Fs);
            end
        end
        
        % calculate RR missed limit
        rr_miss = round(1.8*rr_mean);
        
        % update indices
        last_qrs_i = peak_i;
        ser_back_i = peak_i;
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

function Result = istwave(Signal, candQrs, lastQrs, Fs)
if candQrs-lastQrs < round(0.36*Fs)
    % length of QRS waves
    len = round(0.05*Fs);
    % max slope of the candidate waveform
    begin1 = max(1,candQrs-len);
    slope1 = max(diff(Signal(begin1:candQrs)));
    % max slope of previous QRS waveform
    begin2 = max(1,lastQrs-len);
    slope2 = max(diff(Signal(begin2:lastQrs)));
    % max slope less then half of that of previous QRS
    Result = slope1 < 0.5*slope2;
else
    Result = false;
end