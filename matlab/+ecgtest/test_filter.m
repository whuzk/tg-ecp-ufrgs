import ecgfilter.*
import ecgfastcode.*
import ecgutilities.*
import ecgmath.*
%close all;

% load ecg
Signal = interpret(EDB.e0116,2);
data = Signal.data - Signal.inival;
dataInt = int16(data);

% apply filters
%{
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
%}

%N = length(Signal);
%t = (0:N-1)/Fs;
%noise = 50*sin(2*pi*60*t)';
%Signal2 = Signal + noise;

%[Y1,Rd,R2,TH1,TH2,RR,delay] = rt_detect_qrs_double(data,Signal.fs,50,4);
%[Y2,Rd2,R22,TH12,TH22,RR2,delay2] = detect_qrs_double(data,Signal.fs,50,4);
%{
[Y1,R,R2,TH1,TH2,RR,delay] = detect_qrs_int(dataInt,Signal.fs);
Y1 = double(Y1);
R = double(R);
R2 = double(R2);
TH1 = double(TH1);
TH2 = double(TH2);
RR = double(RR);
%}

%figure, plot(Signal);
%figure, plot([Y1 Y2]);
%norm(Y1-Y2)

%plot_signal_qrs(Y1,Rd,R2,TH1,TH2,RR);
%Rd = Rd - floor(delay);
%plot_signal_r(data, R);
%[A,B,R] = ecgutilities.merge_qrs(Signal.qrs, Rd, Signal.fs);
%ecgutilities.plot_signal_rcomp(data,A,B,R);
%ecgmath.compute_statistics(A,B)