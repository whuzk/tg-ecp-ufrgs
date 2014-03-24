function [QRSi,QRS2,THRs,THRn,THRs2,THRn2] = tompkins_original(SignalB,SignalI,Fs)

%% global initializations
N = length(SignalI);        % length of MWI signal
Ttrain = 2*Fs;              % length of training period
Tref = round(0.20*Fs);      % length of refractory period
TQtol = round(0.10*Fs);     % tolerance for comparing QRS waves
TTtol = round(0.36*Fs);     % tolerance for T wave identification
TPtol = round(0.05*Fs);     % tolerance for finding peaks in BP signal

%% initializations for phase 1
THRs = zeros(N,1);      % buffer for the Signal Threshold history
THRn = zeros(N,1);      % buffer for the Noise Threshold history
THRs2 = zeros(N,1);     % buffer for the Signal Threshold history
THRn2 = zeros(N,1);     % buffer for the Noise Threshold history
THR_SIG1 = 0;           % Signal threshold in MWI
SIG_LEV1 = 0;           % Signal level in MWI
NOISE_LEV1 = 0;         % Noise level in MWI
THR_SIG2 = 0;           % Signal threshold in Bandpass
SIG_LEV2 = 0;           % Signal level in Bandpass
NOISE_LEV2 = 0;         % Noise level in Bandpass

%% algorithm - phase 1
for i = 2:Ttrain
    
    % check if current sample is a peak
    if ispeak(SignalI, i)
        
        % get current sample
        a = SignalI(i);
        
        % find peak in the bandpass signal
        y = findmax(SignalB, i, TPtol);
        
        % check if peak is from qrs
        if a >= THR_SIG1 && y >= THR_SIG2
            % adjust signal levels
            SIG_LEV1 = 0.125*a + 0.875*SIG_LEV1;
            SIG_LEV2 = 0.125*y + 0.875*SIG_LEV2;
        else
            % adjust noise levels
            NOISE_LEV1 = 0.125*a + 0.875*NOISE_LEV1;
            NOISE_LEV2 = 0.125*y + 0.875*NOISE_LEV2;
        end
    end
    
    % update thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    
    % update threshold history
    THRs(i) = THR_SIG1;
    THRn(i) = 0.5*THR_SIG1;
    THRs2(i) = THR_SIG2;
    THRn2(i) = 0.5*THR_SIG2;
end

%% initializations for phase 2
QRSi = zeros(N,1);      % buffer for the R wave indices
act_search = false;     % flag to indicate the active search
ref_count = 0;          % counter for the refractory period
qrs_count = 0;          % count of QRS complex in main QRS buffer
last_qrs_i = 0;         % index of last detected QRS complex
cur_i = 0;              % index of current candidate peak
cur_a = 0;              % amplitude of current candidate peak in MWI
cur_y = 0;              % amplitude of current candidate peak in BP
new_rr = 0;             % new RR interval after QRS detection

%% algorithm - phase 2
i = Ttrain+1;
while new_rr == 0 && i < length(SignalI)
    
    % get current sample
    a = SignalI(i);
    
    % check if the active search is on
    if act_search && a < 0.5*cur_a
        
        % check if candidate peak is from qrs
        if cur_y >= THR_SIG2
            % increment qrs count
            qrs_count = qrs_count + 1;
            % save index of MWI
            QRSi(qrs_count) = cur_i;
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
    
    elseif ispeak(SignalI, i)
        
        % find peak in the bandpass signal
        y = findmax(SignalB, i, TPtol);
        
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
    
    % update threshold history
    THRs(i) = THR_SIG1;
    THRn(i) = 0.5*THR_SIG1;
    THRs2(i) = THR_SIG2;
    THRn2(i) = 0.5*THR_SIG2;
    
    i = i + 1;
end

%% initializations for phase 3
QRS2 = zeros(N,1);              % buffer for QRS detected in searchback
RR_int = NaN(1,8);              % buffer for previous RR intervals
sel_RR_int = NaN(1,8);          % buffer for selected RR intervals
qrs_count2 = 0;                 % current position in second QRS buffer
rr_i = 1;                       % current position in main RR buffer
sel_rr_i = 0;                   % current position in second RR buffer
RR_int(rr_i) = new_rr;          % initialize main RR buffer
rr_mean = nanmean(RR_int);      % running avereage of RR intervals
rr_miss = round(1.66*rr_mean);  % interval for qrs to be assumed missed
qrs_updated = false;            % flag to indicate detection of QRS
ser_back_i = 0;                 % index of searchback starting point

%% algorithm - phase 3
for i = i:length(SignalI)-1
    
    % get current sample
    a = SignalI(i);
    
    % check if we are in the refractory period
    if ref_count > 0
        
        ref_count = ref_count - 1;
        
    elseif act_search && a < 0.5*cur_a
    
        % check if candidate peak is from qrs
        if cur_y >= THR_SIG2
            % skip when a T wave is detected
            if ~istwave(SignalI, cur_i, last_qrs_i, TTtol, TQtol)
                % increment qrs count
                qrs_count = qrs_count + 1;
                % save index of MWI
                QRSi(qrs_count) = cur_i;
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
        
    elseif ispeak(SignalI, i)
    
        % find peak in the bandpass signal
        y = findmax(SignalB, i, TPtol);
        
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
        [new_a, new_i] = findmax(SignalI, center, round(0.5*rr_mean));
        
        % find peak in the bandpass signal
        new_y = findmax(SignalB, new_i, TPtol);

        % check if candidate peak is from qrs
        if new_a > 0.5*THR_SIG1 && new_y >= 0.5*THR_SIG2
            % skip when a T wave is detected
            if ~istwave(SignalI, new_i, last_qrs_i, TTtol, TQtol)
                % update last QRS peak and RR interval
                new_rr = new_i - last_qrs_i;
                last_qrs_i = new_i;
                qrs_updated = true;
                % increment qrs count
                qrs_count2 = qrs_count2 + 1;
                % save index of MWI
                QRS2(qrs_count2) = new_i;
                % activate refractory period
                ref_count = max(0,Tref-(i-new_i));
                % adjust signal levels
                SIG_LEV1 = 0.25*new_a + 0.75*SIG_LEV1;
                SIG_LEV2 = 0.25*new_y + 0.75*SIG_LEV2;
            end
        end
    end
    
    % check whether a new QRS was detected
    if qrs_updated
        
        % update RR intervals and calculate running average
        rr_i = mod(rr_i,length(RR_int))+1;
        RR_int(rr_i) = new_rr;
        rr_mean = nanmean(RR_int);
        
        % check if new beat is regular
        if new_rr >= 0.92*rr_mean && new_rr <= 1.16*rr_mean
            % update buffer of selected RR intervals
            sel_rr_i = mod(sel_rr_i,length(sel_RR_int))+1;
            sel_RR_int(sel_rr_i) = new_rr;
        end
        
        % check if heart rate is regular
        irregular = RR_int < 0.92*rr_mean | RR_int > 1.16*rr_mean;
        if ~isempty(find(irregular, 1))
            % substitute the average
            rr_mean = nanmean(sel_RR_int);
            % adjust levels
            SIG_LEV1 = 0.5*SIG_LEV1;
            SIG_LEV2 = 0.5*SIG_LEV2;
            NOISE_LEV1 = 0.5*NOISE_LEV1;
            NOISE_LEV2 = 0.5*NOISE_LEV2;
        end
        % calculate RR miss limit
        rr_miss = round(1.66*rr_mean);
        
        % enable searchback
        ser_back_i = last_qrs_i;
        
        % reset detection flag
        qrs_updated = false;
    end
    
    % update thresholds
    THR_SIG1 = NOISE_LEV1 + 0.25*abs(SIG_LEV1 - NOISE_LEV1);
    THR_SIG2 = NOISE_LEV2 + 0.25*abs(SIG_LEV2 - NOISE_LEV2);
    
    % update threshold history
    THRs(i) = THR_SIG1;
    THRn(i) = 0.5*THR_SIG1;
    THRs2(i) = THR_SIG2;
    THRn2(i) = 0.5*THR_SIG2;
end

% trim the output vectors
QRSi(qrs_count+1:end) = [];
QRS2(qrs_count2+1:end) = [];


function Result = ispeak(Signal, i)
% get current and neighbour samples
a = Signal(i-1);
b = Signal(i);
c = Signal(i+1);
Result = a-b < 0 && b-c >= 0;

function [y,x] = findmax(Signal,i,l)
N = length(Signal);
begin = max(1,i-l);
end_i = min(N,i+l);
[y,x] = max(Signal(begin:end_i));
x = x + begin - 1;

function Result = istwave(Signal, candQrs, lastQrs, tolerance, qrsLen)
% check if QRS candidate occurs near the previous QRS
if candQrs-lastQrs <= tolerance
    % max slope of the candidate waveform
    begin1 = max(1,candQrs-qrsLen);
    slope1 = max(diff(Signal(begin1:candQrs)));
    % max slope of previous QRS waveform
    begin2 = max(1,lastQrs-qrsLen);
    slope2 = max(diff(Signal(begin2:lastQrs)));
    % slope less then 0.5 of previous QRS
    Result = slope1 < 0.5*slope2;
else
    Result = false;
end