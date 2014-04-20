[Fs,Bp,~,~,Data] = ecgutilities.interpret(EDB.e0116);
Signal = Data(:,2);
%Signal = linspace(0,1000,1000000);
%Signal = linspace(1000,0,1000000);

total = 250;
Smax = zeros(total,1);
Cmax = zeros(total,1);
Smean = zeros(total,1);
Cmean = zeros(total,1);
for i = 1:total
    [~,mem,cmp] = ecgmath.running_max(Signal,i);
    Smax(i) = max(mem);
    Cmax(i) = max(cmp);
    Smean(i) = mean(mem);
    Cmean(i) = mean(cmp);
    i
end
figure, plot([Smax Cmax]);
figure, plot([Smean Cmean]);