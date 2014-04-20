import ecgfilter.*
import ecgfastcode.*
import ecgutilities.*
import ecgmath.*
%close all;

% load ecg
[Fs,Bp,~,Leads] = interpret(EDB.e0103);
Signal = Leads{1}.data;

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

Ws = floor(0.05*Fs)*2+1;
h1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
h2 = ones(1,Ws)/Ws;

X = abs(wconv1(Signal,h1,'same'));
G = running_max(X,2*Fs);
%G = wconv1(G,ones(1,1000)/1000,'same');
Q = (2^8-1)*min(1,X./G);
I = wconv1(mobd(Q,3),h2,'same');

figure, plot(Signal);
figure, plot([X G]);
figure, plot(Q);
figure, plot(I);

%{
[Rdetected,R2,TH1,TH2,RR] = ecgfastcode.c_detect_qrs(W,Fs);
ecgutilities.plot_signal_qrs(Signal,W,Rdetected,R2,TH1,TH2,RR);
[A,B,R] = ecgutilities.merge_rpeaks(Bp, Rdetected, Fs);
ecgutilities.plot_signal_rcomp(Signal,A,B,R);
ecgmath.compute_statistics(A,B)
%}