function Result = fiducial_marks(beat,lead,Rpeak,Fs)

%
N = length(beat);
S = floor(0.05*Fs);
L1 = floor(0.08*Fs);
L2 = floor(0.1*Fs);
L3 = floor(0.3*Fs);

% transform
Temp = ecgmmo.mmo_derivative_flat(beat,S);

% decide which signals to use
if ismember(lead,{'MLIII'})
    TempL = -del2(Temp);
    TempR = -Temp;
    TempP = Temp;
    TempT = -Temp;
else
    TempL = del2(Temp);
    TempR = Temp;
    TempP = Temp;
    TempT = Temp;
end

% first half
Ron = search_first_mark(TempR,Rpeak-1,-1,1);
k = search_best_mark(TempL,Ron-1,-1,max(1,Ron-L1));
Ppeak = search_best_mark(-TempP,k-5,-1,max(1,k-L3));
Pon = search_best_mark(TempP,Ppeak-1,-1,max(1,Ppeak-L1));
Poff = search_best_mark(TempP,Ppeak+1,1,min(k,Ppeak+L1));

% second half
Roff = search_first_mark(TempR,Rpeak+1,1,N);
k = search_best_mark(TempL,Roff+1,1,min(N,Roff+L1));
Tpeak = search_best_mark(-TempT,k+5,1,min(N,k+L3));
Ton = search_best_mark(TempT,Tpeak-1,-1,max(k,Tpeak-L2));
Toff = search_best_mark(TempT,Tpeak+1,1,min(N,Tpeak+L2));

% save result
Result.P = [Pon Ppeak Poff];
Result.R = [Ron Rpeak Roff];
Result.T = [Ton Tpeak Toff];
%{
ecgutilities.plot_fiducial_marks(beat,Result);
ecgutilities.plot_fiducial_marks(Temp,Result);
ecgutilities.plot_fiducial_marks(TempL,Result);
%}

function pos = search_first_mark(data,istart,inc,iend)
if abs(iend-istart)+1 >= 3
    win = data(istart:inc:iend);
    [~,x] = findpeaks(win,'NPeaks',1);
    if isempty(x)
        [~,x] = max(win);
    end
    pos = istart + inc*(x-1);
else
    pos = iend;
end

function pos = search_best_mark(data,istart,inc,iend)
if abs(iend-istart)+1 >= 3
    win = data(istart:inc:iend);
    [y,x] = findpeaks(win);
    if isempty(x)
        [~,x] = max(win);
    else
        [~,i] = max(y);
        x = x(i);
    end
    pos = istart + inc*(x-1);
else
    pos = iend;
end