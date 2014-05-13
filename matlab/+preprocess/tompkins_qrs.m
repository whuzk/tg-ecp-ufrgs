function [qrs1,qrs2,th1,th2,rr1,rr2] = tompkins_qrs(sigI, Fs)

%% initializations
N = length(sigI);               % length of the signal
th1 = zeros(N,1);               % buffer for the threshold history
th2 = zeros(N,1);               % buffer for the secondary threshold history
rr1 = zeros(N,1);               % buffer for the RR mean history
rr2 = zeros(N,1);               % buffer for the secondary RR mean history
qrs1 = zeros(N,1);              % buffer for the QRS locations
qrs2 = zeros(N,1);              % buffer for QRS detected in searchback

training_period = 2*Fs;         % length of training period
threshold = 0;                  % signal threshold
sig_level = 0;                  % signal level
noise_level = 0;                % noise level
est_ratio = 0.25;               % ratio of signal/noise level estimation
qrs_count = 0;                  % count of QRS complex in main QRS buffer
qrs_count2 = 0;                 % current position in second QRS buffer
searchback_idx = 0;             % index of searchback starting point
last_peak_idx = 0;              % index of the last peak detected
last_peak_amp = 0;              % amplitude of the last peak detected
last_qrs_idx = 0;               % index of last detected QRS complex
rr_mean = 0;                    % running average of RR intervals
rr_mean2 = 0;                   % secondary running average of RR intervals
rr_miss = 0;                    % interval for qrs to be assumed as missed
signal_rising = false;          % flag to indicate a rise in the signal

%% algorithm
for i = 1:N
    
    % get current sample
    a = sigI(i);
    
    % initialize iteration flags
    peak_detected = false;
    qrs_detected = false;
    
    %% peak detection
    if a > last_peak_amp
        
        % signalize positive slope
        signal_rising = true;
        
        % update peak info
        last_peak_amp = a;
        last_peak_idx = i;
        
    elseif ~signal_rising
        
        % update current amplitude
        last_peak_amp = a;
        
    elseif a < 0.5*last_peak_amp
        
        % signalize peak detection
        peak_amp = last_peak_amp;
        peak_idx = last_peak_idx;
        peak_detected = true;
        
        % reset state
        signal_rising = false;
        last_peak_amp = a;
        
    end
    
    %% qrs detection and level estimation
    if peak_detected
        if i <= training_period
            if peak_amp >= threshold
                sig_level = max(sig_level,peak_amp);
                last_qrs_idx = peak_idx;
            else
                noise_level = max(noise_level,peak_amp);
            end
        elseif peak_amp >= threshold
            sig_level = estimate(sig_level,est_ratio,peak_amp);
            qrs_detected = ~istwave(sigI, peak_idx, last_qrs_idx, rr_mean2, Fs);
        else
            noise_level = estimate(noise_level,est_ratio,peak_amp);
        end
    end
    
    %% search back
    if qrs_detected || signal_rising || searchback_idx == 0
        
        % do nothing if:
        %   1. a qrs complex was detected in the current iteration
        %   2. a new qrs might be detected in the following iterations
        %   3. the search back procedure is not yet enabled
        
    elseif i-searchback_idx >= rr_miss
        
        % search back and locate the max in this interval
        interval = floor(rr_mean2);
        peak_idx = findmax(sigI, i-interval+1, interval);
        peak_amp = sigI(peak_idx);
        
        % check if candidate peak is from qrs
        if (peak_amp < threshold) && (peak_amp >= 0.5*threshold) && ...
                ~istwave(sigI, peak_idx, last_qrs_idx, rr_mean2, Fs)
            % increment qrs count
            qrs_count2 = qrs_count2 + 1;
            % save index of integrated signal
            qrs2(qrs_count2) = peak_idx;
            % signalize qrs detection
            qrs_detected = true;
            % adjust signal levels
            sig_level = estimate(sig_level,0.25,peak_amp);
        else
            % reduce levels by half
            sig_level = 0.5*sig_level;
            noise_level = 0.5*noise_level;
            % postpone searchback
            searchback_idx = searchback_idx + interval;
        end
    end
    
    
    %% update QRS info and RR-interval
    if qrs_detected
        
        % push new QRS location to output buffer
        qrs_count = qrs_count + 1;
        qrs1(qrs_count) = peak_idx;
        
        % update RR intervals and limits
        new_rr = peak_idx - last_qrs_idx;
        [rr_mean,rr_mean2] = update_rr(new_rr);
        rr_miss = round(1.66*rr_mean2);
        
        % update indices
        last_qrs_idx = peak_idx;
        searchback_idx = peak_idx;
    end
    
    %% updates
    threshold = noise_level + 0.25*abs(sig_level - noise_level);
    
    % decrease estimation ratio for times beyond the training period
    if i == training_period
        est_ratio = est_ratio/2;
    end

    % history
    th1(i) = threshold;
    th2(i) = 0.5*threshold;
    rr1(i) = rr_mean;
    rr2(i) = rr_mean2;
end

% trim the output vectors
qrs1(qrs_count+1:end) = [];
qrs2(qrs_count2+1:end) = [];


function Result = estimate(x,r,v)
Result = (1-r)*x + r*v;

function Result = findmax(signal, begin, interval)
[~,x] = max(signal(begin:begin+interval-1));
Result = x + begin - 1;

function Result = maxdiff(signal, begin, interval)
Result = max(diff(signal(begin:begin+interval-1)));

function Result = istwave(sigI, candQrs, lastQrs, rrMean, Fs)
rr = candQrs - lastQrs;
if rr <= floor(0.2*Fs)
    Result = true;
elseif rr > floor(0.5*rrMean)
    Result = false;
else
    % half the length of a QRS
    len = round(0.1*Fs);
    % find max slopes
    slope1 = maxdiff(sigI,candQrs-len+1,len);
    slope2 = maxdiff(sigI,lastQrs-len+1,len);
    % check condition for T wave
    Result = slope1 < 0.5*slope2;
end

function [rr_mean,rr_mean2] = update_rr(new_rr)
persistent init rr_last rr_last2 rr_idx rr_idx2 rr1 rr2 rr_low rr_high
if isempty(init)
    rr_idx = 1;
    rr_idx2 = 1;
    rr_last(1:8) = new_rr;
    rr_last2(1:8) = new_rr;
    rr_low = 0.92*new_rr;
    rr_high = 1.16*new_rr;
    rr1 = new_rr;
    rr2 = new_rr;
    init = true;
else
    rr_idx = mod(rr_idx,8) + 1;
    rr1 = rr1 + 0.125*(new_rr - rr_last(rr_idx));
    rr_last(rr_idx) = new_rr;
    if rr_low <= new_rr && new_rr <= rr_high
        rr_idx2 = mod(rr_idx2,8) + 1;
        rr2 = rr2 + 0.125*(new_rr - rr_last2(rr_idx2));
        rr_last2(rr_idx2) = new_rr;
        rr_low = 0.92*rr2;
        rr_high = 1.16*rr2;
    end
end
rr_mean = rr1;
rr_mean2 = rr2;