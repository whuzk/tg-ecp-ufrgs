import ecgfastcode.*
import ecgutilities.*
%close all;

% load ecg
Signal = interpret(EDB.e0103,2);
data = Signal.data - Signal.inival;

% detect qrs
[Rd0,Rint,delay0] = prod_detect_qrs_double(data,Signal.fs,50,3);
[Y1,Rd1,R21,TH11,TH21,RR1,delay1] = rt_detect_qrs_double(data,Signal.fs,50,3);
[Y2,Rd2,R22,TH12,TH22,RR2,delay2] = detect_qrs_double(data,Signal.fs,50,3);

%{
dataInt = int16(data);
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

%plot_signal_qrs(Y1,Rd1,R21,TH11,TH21,RR1);
%Rd = Rd - floor(delay);
%plot_signal_r(data, R);
%[A,B,R] = merge_qrs(Signal.qrs, Rd2, Signal.fs);
%plot_signal_rcomp(data,A,B,R);
%ecgmath.compute_statistics(A,B)