Signal = Database.e0112.Signals{1}.Data;
Signal = Signal - Signal(1);
Fs = 250;

% butterworth
fc = [5 15]*2/Fs;
[b,a] = butter(2,fc);
d = 11;

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];
h_l = filter(b_l, a_l, [1 zeros(1,12)]);

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];
h_h = filter(b_h, a_h, [1 zeros(1,32)]);

% derivative filter (2T delay)
h_d = 0.1*[-1 -2 0 2 1];

% integration filter (15T delay)
Ws = round(0.15*Fs);
h_i = ones(1,Ws)/Ws;

% apply filters and update delays
SignalH1 = filter(b, a, Signal);
SignalH1 = [zeros(d,1); SignalH1(1:end-d)];
SignalD1 = filter(h_d, 1, SignalH1);
SignalI1 = filter(h_i, 1, SignalD1.^2);

SignalL = filter(b_l, a_l, Signal);
SignalH2 = filter(b_h, a_h, SignalL);
SignalD2 = filter(h_d, 1, SignalH2);
SignalI2 = filter(h_i, 1, SignalD2.^2);

figure, plot([SignalI1 SignalI2]);