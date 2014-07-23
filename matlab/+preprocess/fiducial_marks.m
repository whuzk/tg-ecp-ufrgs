function Result = fiducial_marks(sigD,sigM,Rpeaks,Fs)
% Funçao para detectar os pontos caracteristicos das ondas de ECG, com
% base na localizaçao dos picos de onda R. Os pontos detectados sao:
%   [Pon Ppk Ron Rpk Rof Tpk Tof]
%   Pon - inicio da onda P
%   Ppk - pico da onda P
%   Ron - inicio da onda R
%   Rpk - pico da onda R
%   Rof - fim da onda R
%   Tpk - pico da onda T
%   Tof - fim da onda T

% initializations
N = length(sigM);
M = 10;%length(Rpeaks);
RR = diff(Rpeaks);
RR1 = [RR(1); RR];
RR2 = [RR; RR(end)];

% limits
L1 = round(0.10*Fs);
L2 = round(0.02*Fs);
L3 = round(0.15*Fs);
L4 = round(0.12*Fs);

% outputs
Result = zeros(7,M);

% algorithm
for i = 1:M
    Rpk = Rpeaks(3000+i-1);
    
    % R-wave
    Rpk = search_peak_abs(sigD,sigM,max(1,Rpk-L1),1,min(N,Rpk+L1),Rpk);
    if sigM(Rpk) > 0
        % inverted
        thresh = abs(sigM(Rpk))/2;
        Ron = search_first_mark(-sigM,Rpk-L2,-1,max(1,Rpk-L1),Rpk);
        Roff = search_first_mark(-sigM,Rpk+L2,1,min(N,Rpk+L4),Rpk);
        Ron = search_mark_below(sigM,Ron-1,-1,max(1,Rpk-L1),Ron,thresh);
        Roff = search_mark_below(sigM,Roff+1,1,min(N,Rpk+L4),Roff,thresh);
    else
        % normal
        thresh = abs(sigM(Rpk))/2;
        Ron = search_first_mark(sigM,Rpk-L2,-1,max(1,Rpk-L1),Rpk);
        Roff = search_first_mark(sigM,Rpk+L2,1,min(N,Rpk+L4),Rpk);
        Ron = search_mark_below(-sigM,Ron-1,-1,max(1,Rpk-L1),Ron,thresh);
        Roff = search_mark_below(-sigM,Roff+1,1,min(N,Rpk+L4),Roff,thresh);
    end
    
    % P-wave
    len = floor(0.375*RR1(i));
    Ppk = search_peak_abs(sigD,sigM,Ron-1,-1,max(1,Ron-len),Ron);
    if sigM(Ppk) > 0
        % inverted
        Pon = search_best_mark(-sigM,Ppk-1,-1,max(1,Ppk-L3),Ppk);
    else
        % normal
        Pon = search_best_mark(sigM,Ppk-1,-1,max(1,Ppk-L3),Ppk);
    end
    
    % T-wave
    len = floor(0.5*RR2(i));
    Tpk = search_peak_abs(sigD,sigM,Roff+1,1,min(N,Roff+len),Roff);
    if sigM(Tpk) > 0
        % inverted
        Toff = search_best_mark(-sigM,Tpk+1,1,min(N,Tpk+L3),Tpk);
    else
        % normal
        Toff = search_best_mark(sigM,Tpk+1,1,min(N,Tpk+L3),Tpk);
    end
    
    % save result
    Result(:,i) = [Pon Ppk Ron Rpk Roff Tpk Toff]';
end


function pos = search_peak_abs(dataD,dataI,istart,inc,iend,default)
idx = [default default];
val = [0 0];
for i = istart+inc:inc:iend-inc
    y = dataD(i);
    if (dataD(i-1) < y && y >= dataD(i+1))
        if y > val(1)
            idx(1) = i;
            val(1) = y;
        end
    elseif (dataD(i-1) > y && y <= dataD(i+1))
        if y < val(2)
            idx(2) = i;
            val(2) = y;
        end
    end
end
if idx(1) < idx(2)
    left = idx(1);
    right = idx(2);
else
    left = idx(2);
    right = idx(1);
end
if right - left <= 1
    pos = left;
else
    [y,x] = findpeaks(abs(dataI(left:right)));
    if ~isempty(x)
        [~,i] = max(y);
        x = x(i);
        pos = left + x - 1;
    else
        pos = left;
    end
end


function pos = search_first_mark(data,istart,inc,iend,default)
if abs(iend - istart - inc) > 0
    win = data(istart:inc:iend);
    [~,x] = findpeaks(win,'NPeaks',1);
    if isempty(x)
        %[~,x] = max(win);
        pos = default;
        return;
    end
    pos = istart + inc*(x-1);
else
    pos = default;
end

function pos = search_mark_below(data,istart,inc,iend,default,thresh)
pos = default;
if abs(iend - istart - inc) > 0
    for i = istart+inc:inc:iend-inc
        y = data(i);
        if data(i-1) < y && y >= data(i+1) && abs(y) < thresh
            pos = i;
            break;
        end
    end
end

function pos = search_best_mark(data,istart,inc,iend,default)
if abs(iend - istart - inc) > 0
    win = data(istart:inc:iend);
    [y,x] = findpeaks(win);
    if ~isempty(y)
        [~,i] = max(y);
        x = x(i);
    else
        [~,x] = max(win);
    end
    pos = istart + inc*(x-1);
else
    pos = default;
end