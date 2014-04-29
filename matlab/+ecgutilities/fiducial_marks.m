function Result = fiducial_marks(mmd,lap,lead,Rpeaks,Fs)

% decide which signals to use
if ismember(lead,{'MLIII'})
    TempL = -lap;
    TempR = -mmd;
    TempP = mmd;
    TempT = -mmd;
else
    TempL = lap;
    TempR = mmd;
    TempP = mmd;
    TempT = mmd;
end

%
N = length(mmd);
M = length(Rpeaks);

%
L1 = floor(0.08*Fs);
L2 = floor(0.1*Fs);
L3 = floor(0.3*Fs);
L4 = floor(0.4*Fs);

P = cell(M,1);
R = cell(M,1);
T = cell(M,1);
for i = 1:M
    Rpeak = Rpeaks(i);
    
    % first half
    Ron = search_first_mark(TempR,Rpeak-1,-1,max(1,Rpeak-L4),0);
    k = search_best_mark(TempL,Ron-1,-1,max(1,Ron-L1));
    Ppeak = search_best_mark(-TempP,k-5,-1,max(1,k-L3));
    Pon = search_best_mark(TempP,Ppeak-1,-1,max(1,Ppeak-L1));
    Poff = search_best_mark(TempP,Ppeak+1,1,min(k,Ppeak+L1));

    % second half
    Roff = search_first_mark(TempR,Rpeak+1,1,min(N,Rpeak+L4),0);
    k = search_best_mark(TempL,Roff+1,1,min(N,Roff+L1));
    Tpeak = search_best_mark(-TempT,k+5,1,min(N,k+L4));
    Ton = search_best_mark(TempT,Tpeak-1,-1,max(k,Tpeak-L2));
    Toff = search_best_mark(TempT,Tpeak+1,1,min(N,Tpeak+L2));
    
    P{i} = [Pon Ppeak Poff];
    R{i} = [Ron Rpeak Roff];
    T{i} = [Ton Tpeak Toff];
end
Result = struct('P',P,'R',R,'T',T);


function pos = search_first_mark(data,istart,inc,iend,thr)
if abs(iend-istart)+1 >= 3
    win = data(istart:inc:iend);
    [~,x] = findpeaks(win,'NPeaks',1,'MinPeakHeight',thr);
    if isempty(x)
        [y,x] = max(win);
        if y < thr
            x = 1;
        end
    end
    pos = istart + inc*(x-1);
else
    pos = istart;
end

function pos = search_best_mark(data,istart,inc,iend)
if abs(iend-istart)+1 >= 3
    win = data(istart:inc:iend);
    [y1,x1] = findpeaks(win);
    if ~isempty(y1)
        [~,i] = max(y1);
        y1 = y1(i);
        x1 = x1(i);
    end
    [y2,x2] = max(win);
    if isempty(y1) || y1 < y2
        x = x2;
    else
        x = x1;
    end
    pos = istart + inc*(x-1);
else
    pos = istart;
end