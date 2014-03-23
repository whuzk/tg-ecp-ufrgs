function QRS = tompkins_production(SignalB, SignalI, DelayI, Fs)

%% initializations
N = length(SignalI);    % length of MWI signal
Tref = round(0.2*Fs);   % length of refractory period
Ttrain = 2*Fs;          % length of training period
Tqrs = round(0.05*Fs);  % half the length of a QRS wave
TTtol = round(0.36*Fs); % tolerance for T wave identification

QRS = zeros(N,1);       % buffer for the R wave indices
RR_int = NaN(1,8);      % buffer for the last RR intervals
sel_RR_int = NaN(1,8);  % buffer for the selected RR intervals

ser_back = false;       % flag to indicate the search back
acthung = false;        % flag to indicate the active search
qrs_updated = false;    % flag to indicate when a new qrs is detected
rr_miss = 0;            % interval beyond which to assume a qrs missed
rr_mean = 0;            % running average of latest RR intervals
rr_count = 0;           % count of selected RR intervals in RR buffer
sel_rr_count = 0;       % count of selected RR intervals in RR buffer
qrs_count = 0;          % count of QRS complex in output QRS buffer
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

%% algorithm - phase 1
for i = DelayI+(1:Ttrain)
    
    % check if current sample is a peak
    if ispeak(SignalI, i)
        
        % get current sample
        a = SignalI(i);
        % find peak in the bandpass signal
        y = findmax(SignalB, i-DelayI, DelayI);
        
        % adjust levels and thresholds for MWI signal
        if a >= THR_SIG1
            SIG_LEV1 = 0.125*a + 0.875*SIG_LEV1;
        elseif a >= THR_NOISE1
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
        end
        THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
        THR_NOISE1 = 0.5*THR_SIG1;
        
        % adjust levels and thresholds for bandpass signal
        if y >= THR_SIG2
            SIG_LEV2 = 0.125*y + 0.875*SIG_LEV2;
        elseif y >= THR_NOISE2
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
        THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
        THR_NOISE2 = 0.5*THR_SIG2;
    end
end

%% algorithm - phase 2
i = DelayI+Ttrain+1;
while rr_count == 0 && i < length(SignalI)
    
    % check if current sample is a peak
    if ispeak(SignalI, i)
        
        % get current sample
        a = SignalI(i);
        % find peak in the bandpass signal
        y = findmax(SignalB, i-DelayI, DelayI);
        
        if acthung && a < 0.5*cur_a
        
            if cur_y >= THR_NOISE2
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of MWI
                QRS(qrs_count) = cur_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-cur_i));
                % update last QRS index and RR interval
                if last_qrs_i > 0
                    new_rr = cur_i - last_qrs_i;
                    rr_count = mod(rr_count,length(RR_int))+1;
                    RR_int(rr_count) = new_rr;
                    rr_mean = nanmean(RR_int);
                end
                last_qrs_i = cur_i;
            end
            % reset flags
            acthung = false;
        
        elseif a >= THR_SIG1
            % check if peak is higher
            if ~acthung || a > cur_a
                cur_a = a;
                cur_i = i;
                cur_y = y;
                acthung = true;
            end
            SIG_LEV1 = 0.125*a + 0.875*SIG_LEV1;
        elseif a >= THR_NOISE1
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
        end
        THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
        THR_NOISE1 = 0.5*THR_SIG1;
        
        % adjust levels and thresholds for bandpass signal
        if y >= THR_SIG2
            SIG_LEV2 = 0.125*y + 0.875*SIG_LEV2;
        elseif y >= THR_NOISE2
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
        THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
        THR_NOISE2 = 0.5*THR_SIG2;
    end
    
    i = i + 1;
end

%% algorithm - phase 3
for i = i:length(SignalI)-1
    
    % get current and neighbour samples
    a = SignalI(i);
    d = i-last_qrs_i;
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif acthung && a < 0.5*cur_a
    
        if cur_y >= THR_NOISE2
            % skip when a T wave is detected
            if ~istwave(SignalI, cur_i, last_qrs_i, TTtol, Tqrs)
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of MWI
                QRS(qrs_count) = cur_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-cur_i));
                % update last QRS index and RR interval
                new_rr = cur_i - last_qrs_i;
                last_qrs_i = cur_i;
                qrs_updated = true;
            end
        end
        % reset flag
        acthung = false;
        
    elseif ispeak(SignalI, i)
    
        % find peak in the bandpass signal
        y = findmax(SignalB, i-DelayI, DelayI);
        
        % find noise and QRS peaks
        if a >= THR_SIG1
            % check if peak is higher
            if ~acthung || a > cur_a
                cur_a = a;
                cur_i = i;
                cur_y = y;
                acthung = true;
            end
            SIG_LEV1 = 0.125*a + 0.875*SIG_LEV1;
        elseif a >= THR_NOISE1 && d < round(0.8*rr_mean)
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
        end
        
        % adjust levels for bandpass signal
        if y >= THR_SIG2
            SIG_LEV2 = 0.125*y + 0.875*SIG_LEV2;
        elseif y >= THR_NOISE2
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
        
    elseif ser_back && (a < THR_NOISE1) && d >= rr_miss
        
        % search back and locate the max in this interval
        center = last_qrs_i + round(rr_mean);
        [new_a, new_i] = findmax(SignalI, center, round(0.5*rr_mean));
        
        % find peak in the bandpass signal
        new_y = findmax(SignalB, new_i-DelayI, DelayI);
        
        % skip if it is a T wave
        if new_y >= THR_NOISE2
            % skip when a T wave is detected
            if ~istwave(SignalI, new_i, last_qrs_i, TTtol, Tqrs)
                % update last QRS peak and RR interval
                new_rr = new_i - last_qrs_i;
                last_qrs_i = new_i;
                qrs_updated = true;
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of MWI
                QRS(qrs_count) = new_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-new_i));
                %when found with the second threshold
                SIG_LEV1 = 0.25*new_a + 0.75*SIG_LEV1;
                % when found with the second threshold
                SIG_LEV2 = 0.25*new_y + 0.75*SIG_LEV2;
            end
        end
        
        % reset flag
        ser_back = false;
        
    elseif d >= 2*rr_miss && mod(d,round(0.2*Fs)) == 0
        
        SIG_LEV1 = 0.5*SIG_LEV1;
        NOISE_LEV1 = 0.5*NOISE_LEV1;
        SIG_LEV2 = 0.5*SIG_LEV2;
        NOISE_LEV2 = 0.5*NOISE_LEV2;
        ser_back = true;
        
    end
    
    % check whether a new QRS was detected and update heart rates
    if qrs_updated
        
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
        end
        
        % calculate RR miss limit
        rr_miss = round(1.66*rr_mean);
        
        % reset flags
        qrs_updated = false;
        ser_back = true;
    end
    
    %% adjust thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_NOISE1 = 0.5*THR_SIG1;

    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    THR_NOISE2 = 0.5*THR_SIG2;
end
% trim the output vector
QRS(qrs_count+1:end) = [];


function Result = ispeak(Signal, i)
% get current and neighbour samples
a = Signal(i-1);
b = Signal(i);
c = Signal(i+1);
Result = a-b < 0 && b-c >= 0;

function [y,x] = findmax(Signal,i,l)
begin = max(1,i-l);
[y,x] = max(Signal(begin:i+l));
x = x + begin - 1;

function Result = istwave(Signal, candQrs, lastQrs, tolerance, qrsLen)
% check if QRS candidate occurs near the previous QRS
if candQrs-lastQrs <= tolerance
    
    % mean slope of the candidate waveform
    begin1 = max(1,candQrs-qrsLen);
    end1 = candQrs;%min(N,candQrs+qrsLen);
    slope1 = max(diff(Signal(begin1:end1)));
    
    % mean slope of previous QRS waveform
    begin2 = max(1,lastQrs-qrsLen);
    end2 = lastQrs;%min(N,lastQrs+qrsLen);
    slope2 = max(diff(Signal(begin2:end2)));
    
    % slope less then 0.5 of previous R
    Result = abs(slope1) < abs(0.5*slope2);
else
    Result = false;
end