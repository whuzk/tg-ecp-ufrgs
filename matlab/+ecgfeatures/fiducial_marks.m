function Result = fiducial_marks(sigD,sigI,Rpeaks,RR,Fs)
% Funçao para detectar os pontos caracteristicos das ondas de ECG, com
% base na localizaçao dos picos de onda R. Os pontos detectados sao:
%   [Pon Ppk]       - inicio e pico da onda P
%   [Ron Rpk Rof]   - inicio, pico e fim da onda R
%   [Tpk Tof]       - pico e fim da onda T
%   [Ipt Jpt]       - pontos isoeletrico e J

% initializations
N = length(sigI);
M = 10;%length(Rpeaks);

% limits
L1 = round(0.10*Fs);
L2 = round(0.04*Fs);
L3 = round(0.20*Fs);
L4 = round(0.02*Fs);

% outputs
P = zeros(M,2);
R = zeros(M,3);
T = zeros(M,2);
IJ = zeros(M,2);

% algorithm
for i = 1:M
    Rpk = Rpeaks(i);
    
    % R-wave
    Rpk = search_rpeak(abs(sigD),abs(sigI),max(1,Rpk-L1),min(N,Rpk+L1),Rpk);
    if sigI(Rpk) > 0
        % inverted
        Ron = search_first_mark(-sigI,Rpk-L4,-1,max(1,Rpk-L1),Rpk);
        Roff = search_first_mark(-sigI,Rpk+L4,1,min(N,Rpk+L1),Rpk);
        Ipt = search_first_mark(sigD,Ron-1,-1,max(1,Ron-L2),Ron);
        Jpt = search_first_mark(-sigD,Roff+1,1,min(N,Roff+L2),Roff);
    else
        % normal
        Ron = search_first_mark(sigI,Rpk-L4,-1,max(1,Rpk-L1),Rpk);
        Roff = search_first_mark(sigI,Rpk+L4,1,min(N,Rpk+L1),Rpk);
        Ipt = search_first_mark(-sigD,Ron-1,-1,max(1,Ron-L2),Ron);
        Jpt = search_first_mark(sigD,Roff+1,1,min(N,Roff+L2),Roff);
    end
    
    % P-wave
    len = floor(0.25*RR(i));
    Ppk = search_best_mark(abs(sigI),Ron-L2,-1,max(1,Ron-len),Ron);
    if sigI(Ppk) > 0
        % inverted
        Pon = search_best_mark(-sigI,Ppk-L2,-1,max(1,Ppk-L3),Ppk);
    else
        % normal
        Pon = search_best_mark(sigI,Ppk-L2,-1,max(1,Ppk-L3),Ppk);
    end
    
    % T-wave
    len = floor(0.35*RR(i));
    Tpk = search_best_mark(abs(sigI),Roff+L2,1,min(N,Roff+len),Roff);
    if sigI(Tpk) > 0
        % inverted
        Toff = search_best_mark(-sigI,Tpk+L2,1,min(N,Tpk+L3),Tpk);
    else
        % normal
        Toff = search_best_mark(sigI,Tpk+L2,1,min(N,Tpk+L3),Tpk);
    end
    
    % save result
    P(i,:) = [Pon Ppk];
    R(i,:) = [Ron Rpk Roff];
    T(i,:) = [Tpk Toff];
    IJ(i,:) = [Ipt Jpt];
end
Result = struct('P',P,'R',R,'T',T,'IJ',IJ);


function pos = search_rpeak(dataD,dataI,istart,iend,default)
idx = [default default];
val = [0 0];
j = 1;
for i = istart+1:iend-1
    y = dataD(i);
    if (dataD(i-1) < y && y >= dataD(i+1))
        if y > val(j)
            idx(j) = i;
            val(j) = y;
            j = mod(j,2) + 1;
        end
    end
end
if idx(1) < idx(2)
    [~,x] = max(dataI(idx(1):idx(2)));
    left = idx(1);
else
    [~,x] = max(dataI(idx(2):idx(1)));
    left = idx(2);
end
pos = left + x - 1;

function pos = search_first_mark(data,istart,inc,iend,default)
if abs(iend - istart - inc) > 0
    win = data(istart:inc:iend);
    [~,x] = findpeaks(win,'NPeaks',1);
    if isempty(x)
        [~,x] = max(win);
    end
    pos = istart + inc*(x-1);
else
    pos = default;
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