import ecgfilter.*
import ecgfastcode.*
import ecgutilities.*
import ecgmath.*
%close all;

% load ecg
%Signal = interpret(EDB.e0103,2);
%data = Signal.data - Signal.inival;
%dataInt = int16(data);

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

%[Rd0,Rint] = prod_detect_qrs_double(data,Signal.fs,50,3);
%[Y1,Rd,R2,TH1,TH2,RR,delay] = rt_detect_qrs_double(data,Signal.fs,50,3);
%[Y2,Rd2,R22,TH12,TH22,RR2,delay2] = detect_qrs_double(data,Signal.fs,50,3);

%{
[Rd0,RR0] = prod_detect_qrs_int(dataInt,Signal.fs,50,3);
Rd0 = double(Rd0);
RR0 = double(RR0);
[Y1,Rd1,R21,TH11,TH21,RR1,delay1] = rt_detect_qrs_int(dataInt,Signal.fs,50,3);
Y1 = double(Y1);
Rd1 = double(Rd1);
R21 = double(R21);
TH11 = double(TH11);
TH21 = double(TH21);
RR1 = double(RR1);
[Y2,Rd2,R22,TH12,TH22,RR2,delay2] = detect_qrs_int(dataInt,Signal.fs,50,3);
Y2 = double(Y2);
Rd2 = double(Rd2);
R22 = double(R22);
TH12 = double(TH12);
TH22 = double(TH22);
RR2 = double(RR2);
%}

%figure, plot(Signal);
%figure, plot([Y1 Y2]);
%norm(Y1-Y2)

%plot_signal_qrs(Y1,Rd,R2,TH1,TH2,RR);
%Rd = Rd - floor(delay);
%plot_signal_r(data, R);
%[A,B,R] = merge_qrs(Signal.qrs, Rd2, Signal.fs);
%plot_signal_rcomp(data,A,B,R);
%ecgmath.compute_statistics(A,B)