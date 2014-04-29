function [MMD,LAP,delay] = sogari_filter2(Signal,Fs,Ms)

Signal = Signal(:);
delay = 0;

% low-pass filtering
[Temp,d] = lpf(Signal,Fs,Ms);
delay = delay + d;

% multiscale morphological derivative
[MMD,d] = mmd(Temp,Fs);
delay = delay + d;

% discrete laplacian
[LAP,d] = lap(MMD);
delay = delay + d;

MMD = [zeros(d,1); MMD(1:end-d)];
%{
figure, plot(Temp1);
figure, plot(Temp2);
%}

function [Result,delay] = lpf(Signal,Fs,Ms)
N = round(Fs/Ms);
h = triangularPulse(-N,N,-N+1:N-1);
a = [1 -2 1];
b = conv(h,a)/N;
Result = filter(b,a,Signal);
delay = N-1;

function [Result,delay] = mmd(Signal,Fs)
s = floor(0.06*Fs);
B = ones(2*s+1,1);
Result = (imerode(Signal,B)+imdilate(Signal,B)-2*Signal)./s;
Result = [zeros(s,1); Result(1:end-s)];
delay = s;

function [Result,delay] = lap(Signal)
h = [1 -2 1]/4;
Result = filter(h,1,Signal);
delay = 1;