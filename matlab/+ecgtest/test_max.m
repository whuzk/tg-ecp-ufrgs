[Fs,Bp,~,Leads] = ecgutilities.interpret(EDB.e0305);
Signal = Leads{1}.data;
%Signal = linspace(0,1000,1000000);
%Signal = linspace(1000,0,1000000);

total = 250;
Smax = zeros(total,1);
Cmax = zeros(total,1);
Smean = zeros(total,1);
Cmean = zeros(total,1);
for i = 1:total
    [~,mem,cmp] = ecgmath.running_max(Signal,i,-Inf);
    Smax(i) = max(mem);
    Cmax(i) = max(cmp);
    Smean(i) = mean(mem);
    Cmean(i) = mean(cmp);
    i
end
figure, plot([Smax Cmax]);
figure, plot([Smean Cmean]);