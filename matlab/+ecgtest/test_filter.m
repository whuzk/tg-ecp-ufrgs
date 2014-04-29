import ecgfilter.*
import ecgmath.*
%close all;

% load ecg
Signal = interpret(EDB.e0103,1);
data = Signal.data - Signal.inival;

% apply filters
Y1 = suppress_noise(Signal,Fs);
Y2 = baseline_filter(Signal,Fs);
Y3 = chu_filter(Signal,Fs);
Y4 = krishnan_filter(Signal,Fs);
Y5 = trahanias_filter(Signal,Fs);
Y6 = sun_filter(Signal,Fs);
[~,Y7] = tompkins_filter(Signal, Fs);

% plot figures
figure, plot(Signal);
figure, plot(Y1);
figure, plot(Y2);
figure, plot(Y3);
figure, plot(Y4);
figure, plot(Y5);
figure, plot(Y6);
figure, plot(Y7);