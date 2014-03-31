import ecgfilter.*
close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(1:10000,2);

% load wavelet filter
wname = deblankl('db3');
[~,fname] = wavemngr('fields',wname,'type','file');
F = feval(fname,wname);
W = F/sum(F);

% apply filters
[SignalF1,SignalI1] = tompkins_preprocess(Signal, Fs);
[SignalF2,SignalI2] = tompkins_preprocess_original(Signal, Fs);
SignalI3 = chen_preprocess(Signal, Fs, W);

% plot figures
figure, plot(Signal);
figure, plot(SignalI1);
figure, plot(SignalI2);
figure, plot(SignalI3);