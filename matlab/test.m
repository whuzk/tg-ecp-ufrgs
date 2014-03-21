%Signal = Database.e0103.Signals{1}.Data;
%Signal = Signal - Signal(1);
%Fs = 250;

tic;
Res1 = resample(Signal, 200, Fs, 10, 5);
toc;
%tic;
%Res2 = resample2(Signal, 200, Fs, 10, 5);
%toc;
tic;
Res3 = resample3(Signal, 200, Fs, 10, 5);
toc;

%norm(Res1 - Res2)
%figure, plot([Res1 Res2]);
norm(Res1 - Res3)
figure, plot([Res1 Res3]);