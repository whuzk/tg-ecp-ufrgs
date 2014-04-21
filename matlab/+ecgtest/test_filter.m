import ecgfilter.*
import ecgfastcode.*
import ecgutilities.*
import ecgmath.*
%close all;

% load ecg
[Fs,Bp,~,Leads] = interpret(EDB.e0611);
Signal = Leads{2}.data;%(1:1000);
Signal = (Signal - Signal(1));
SignalInt = int16(Signal);

% apply filters
%{
Y1 = suppress_noise(Signal,Fs);
Y2 = baseline_filter(Signal,Fs);
Y3 = chu_filter(Signal,Fs);
Y4 = krishnan_filter(Signal,Fs);
Y5 = trahanias_filter(Signal,Fs);
Y6 = sun_filter(Signal,Fs);
[~,Y7] = tompkins_filter(Signal, Fs);
Y8 = c_filter_double(Signal,Fs);

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
%}

[Filt,delay] = sogari_filter(Signal,Fs,4);
delay = floor(delay);
%{
%figure, plot(Signal);
[Filt2,delay2] = filter_double(Signal,Fs,3);
%Filt2 = double(filter_int(SignalInt,Fs));
figure, plot(abs(Filt-Filt2));
figure, plot([Filt Filt2]);
norm(Filt-Filt2)
delay-delay2
%}
%
[Rdetected,R2,TH1,TH2,RR] = ecgfastcode.detect_qrs(Filt,Fs);
ecgutilities.plot_signal_r(Signal,Rdetected-delay);
ecgutilities.plot_signal_qrs(Filt,Rdetected,R2,TH1,TH2,RR);
[A,B,R] = ecgutilities.merge_rpeaks(Bp, Rdetected-delay, Fs);
ecgutilities.plot_signal_rcomp(Signal,A,B,R);
ecgmath.compute_statistics(A,B)
%