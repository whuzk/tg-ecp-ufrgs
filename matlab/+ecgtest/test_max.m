Signal = ecgutilities.interpret(EDB.e0116,2);
data = Signal.data - Signal.inival;
%data = data(1:10000);
%data = resample(data, 15360, Signal.fs);
%Fs = 15360;
Fs = Signal.fs;

N = round(Fs/60);
h = triangularPulse(-N,N,-N+1:N-1);
a = [1 -2 1];
b = conv(h,a)/N;
Y = filter(b,a,data);

s = floor(0.06*Fs);
G1 = ecgmath.running_max(Y,2*s+1,-Inf);
G2 = -ecgmath.running_max(-Y,2*s+1,-Inf);
temp1 = (G1+G2-2*[zeros(s,1); Y(1:end-s)])./s;

B = ones(2*s+1,1);
temp2 = (imerode(Y,B)+imdilate(Y,B)-2*Y)./s;
temp2 = [zeros(s,1); temp2(1:end-s)];
figure, plot([temp1 temp2]);

%data = linspace(0,1000,1000000);
%data = linspace(1000,0,1000000);
%{
w = floor(0.06*Fs)*2+1;
[G,mem,cmp] = ecgmath.running_max(-Y,w,-Inf);
max(mem)
max(cmp)
mean(mem)
mean(cmp)
figure, plot(data);
figure, plot([Y -G]);
%}
%{
total = 2*15360;
Smax = zeros(total,1);
Cmax = zeros(total,1);
Smean = zeros(total,1);
Cmean = zeros(total,1);
for i = 15000:total
    [~,mem,cmp] = ecgmath.running_max(Signal,i,-Inf);
    Smax(i) = max(mem);
    Cmax(i) = max(cmp);
    Smean(i) = mean(mem);
    Cmean(i) = mean(cmp);
    i
end
figure, plot([Smax Cmax]);
figure, plot([Smean Cmean]);
%}