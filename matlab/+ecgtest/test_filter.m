import ecgfilter.*
%close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);

% remove noise
SignalD = suppress_noise(Signal,Fs);
%SignalD = wavelet_filter(Signal,Fs,'bior3.5');

% apply filters
[SignalF,SignalI] = tompkins_filter(Signal, Fs);
%R = ecgfeatures.tompkins_production(SignalF, SignalI, Fs);
%SignalD = myfilter(Signal,R);

% plot figures
figure, plot(Signal);
figure, plot(SignalD);
figure, plot(SignalI);