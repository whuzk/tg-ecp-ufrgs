[Fs,Bp,~,Leads] = ecgutilities.interpret(EDB.e0305);
Signal = Leads{2}.data(1:100000);
Signal = (Signal - Signal(1));
Signal = resample(Signal, 15360, Fs);
Fs = 15360;
N = round(Fs/60);
h1 = N*triangularPulse(-N,N,-N+1:N-1);
h2 = [1 -2 1];
h3 = conv(h1,h2)/4;
Y = wconv1(Signal,h3,'same');
%Signal = linspace(0,1000,1000000);
%Signal = linspace(1000,0,1000000);

[G,mem,cmp] = ecgmath.running_max(Y,2*Fs,-Inf);
max(mem)
max(cmp)
mean(mem)
mean(cmp)
figure, plot(Signal);
figure, plot([Y G]);

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