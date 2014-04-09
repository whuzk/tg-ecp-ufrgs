import ecgfilter.*
close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0116);
Signal = Data(180001:190000,2);

% de-noise
SignalD = wden(Signal,'sqtwolog','h','sln',3,'db5');
%SignalD = cyclicrtdenoise3(Signal,128,'db5');

% apply filters
[SignalF1,SignalI1] = tompkins_preprocess(Signal, Fs);
[SignalF2,SignalI2] = tompkins_preprocess(SignalD, Fs);

% plot figures
figure, plot(Signal);
figure, plot(SignalD);
figure, plot(SignalI1);
figure, plot(SignalI2);