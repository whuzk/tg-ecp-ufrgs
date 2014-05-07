function Result = fiducial_marks(sigD,thrD,sigI,lead,Rpeaks,RR,Fs)
% Funçao para detectar os pontos caracteristicos das ondas de ECG, com
% base na localizaçao dos picos de onda R. Os pontos detectados sao:
%   [Pon Ppk]   - inicio e pico da onda P
%   [Ron Rof]   - inicio e fim da onda R
%   [Tpk Tof]   - pico e fim da onda T
%   [Ipt Jpt]   - pontos isoeletrico e J
global pwavepolarity rwavepolarity twavepolarity

% decide which signals to use for R and T waves
TempP = pwavepolarity(lead) * sigI;
TempR = rwavepolarity(lead) * sigI;
TempT = twavepolarity(lead) * sigI;

% initializations
N = length(sigI);
M = length(Rpeaks);

% limits
L0 = round(0.04*Fs);
L1 = round(0.10*Fs);
L2 = round(0.20*Fs);
T1 = round(0.12*Fs);
T2 = round(0.04*Fs);
T3 = round(0.02*Fs);
T4 = round(0.12*Fs);

% outputs
P = zeros(M,2);
R = zeros(M,2);
T = zeros(M,2);
IJ = zeros(M,2);

% algorithm
for i = 1:M
    % get new R peak
    Rpk = Rpeaks(i);
    
    % calculate RR related limits
    L3 = round(0.25*RR(i));
    L4 = round(0.35*RR(i));
    
    % first half of the beat
    Ron = search_first_mark(TempR,Rpk-5,-1,max(1,Rpk-L1),Rpk);
    Ppk = search_best_mark(-TempP,Ron-L0,-1,max(1,Ron-L3),Ron);
    Pon = search_best_mark(TempP,Ppk-L0,-1,max(1,Ppk-L2),Ppk);
    
    % second half of the beat
    Roff = search_first_mark(TempR,Rpk+5,1,min(N,Rpk+L1),NaN);
    Tpk = search_best_mark(-TempT,Roff+L0,1,min(N,Roff+L4),Roff);
    Toff = search_best_mark(TempT,Tpk+L0,1,min(N,Tpk+L2),Tpk);
    
    % isoeletric and J points
    Ipt = search_mark(sigD,Rpk-T2,-1,max(1,Rpk-T1),thrD,Ron);
    Jpt = search_mark(sigD,Rpk+T3,1,min(N,Rpk+T4),thrD,Roff);
    
    % save result
    P(i,:) = [Pon Ppk];
    R(i,:) = [Ron Roff];
    T(i,:) = [Tpk Toff];
    IJ(i,:) = [Ipt Jpt];
end
Result = struct('P',P,'R',R,'T',T,'IJ',IJ);


function pos = search_first_mark(data,istart,inc,iend,default)
if inc*(iend-istart) > 1
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
if inc*(iend-istart) > 1
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

function pos = search_mark(data,istart,inc,iend,thr,default)
if inc*(iend-istart) >= 0
    win = data(istart:inc:iend);
    x = find(abs(win) < thr, 1, 'first');
    if isempty(x)
        [~,x] = min(abs(win));
    end
    pos = istart + inc*(x-1);
else
    pos = default;
end