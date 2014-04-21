import ecgfilter.*
import ecgfastcode.*
import ecgutilities.*
import ecgmath.*
%close all;

% load ecg
[Fs,Bp,~,Leads] = interpret(EDB.e0119);
Signal = Leads{2}.data(1:1000);
Signal = (Signal - Signal(1));
SignalInt = int16(Signal);

% apply filters
%
Y1 = suppress_noise(Signal,Fs);
Y2 = baseline_filter(Signal,Fs);
Y3 = chu_filter(Signal,Fs);
Y4 = krishnan_filter(Signal,Fs);
Y5 = trahanias_filter(Signal,Fs);
Y6 = sun_filter(Signal,Fs);
[~,Y7] = tompkins_filter(Signal, Fs);
Y8 = filter_double(Signal,Fs,3);
Y9 = double(filter_int(SignalInt,Fs,3));

% plot figures
figure, plot(Signal);
figure, plot(Y1);
figure, plot(Y2);
figure, plot(Y3);
figure, plot(Y4);
figure, plot(Y5);
figure, plot(Y6);
figure, plot(Y7);
figure, plot(Y8);
figure, plot(Y9);
%