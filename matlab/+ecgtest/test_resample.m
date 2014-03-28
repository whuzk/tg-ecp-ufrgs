import ecgfilter.*

%% initialization
Signal = Database.e0103.Signals{1}.Data;
Signal = Signal - Signal(1);

%% filter design
[p,q] = rat(NewFs/Fs, 1E-12);
pqmax = max(p,q);
fc = 1/2/pqmax;
L = 2*10*pqmax+1;
f_r = [0 2*fc 2*fc 1];
a_r = [1 1 0 0];
h_r = firls(L-1, f_r, a_r);
h_r = p*h_r.*kaiser(L,5)';

%% resampling
tic;
Res1 = resample(Signal, NewFs, Fs, h_r);
toc;
tic;
Res2 = resample2(Signal, NewFs, Fs, h_r);
toc;
tic;
Res3 = resample3(Signal, NewFs, Fs, h_r);
toc;

%% verification
norm(Res1 - Res2)
figure, plot([Res1 Res2]);
norm(Res1 - Res3)
figure, plot([Res1 Res3]);

%% inverse resampling
tic;
InvRes1 = resample(Res1, Fs, NewFs, h_r*q/p);
toc;
tic;
InvRes2 = resample2(Res2, Fs, NewFs, h_r*q/p);
toc;
tic;
InvRes3 = resample3(Res3, Fs, NewFs, h_r*q/p);
toc;

%% verification
norm(Signal - InvRes1)
figure, plot([Signal InvRes1]);
norm(Signal - InvRes2)
figure, plot([Signal InvRes2]);
norm(Signal - InvRes3)
figure, plot([Signal InvRes3]);
