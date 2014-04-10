import ecgfilter.*
close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);

% remove noise
SignalD = wavfilter(Signal,round(log2(Fs)),'db5');

% apply filters
[SignalF1,SignalI1] = tompkins_preprocess(Signal, Fs);
[SignalF2,SignalI2] = tompkins_preprocess(SignalD, Fs);

% plot figures
figure, plot(Signal);
figure, plot(SignalD);
figure, plot(SignalI1);
figure, plot(SignalI2);