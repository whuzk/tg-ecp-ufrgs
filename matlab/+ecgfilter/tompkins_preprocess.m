function [Result,Delay] = tompkins_preprocess(Signal, Fs)

% low-pass (5T delay)
A1 = [1 -2 1];
B1([1 7 13]) = 1/36*[1 -2 1];
T1 = tf(B1,A1);

% high-pass (16T delay)
A2 = [1 -1];
B2([1 17 18 33]) = [-1/32 1 -1 1/32];
T2 = tf(B2,A2);

% derivative (2T delay)
A3 = 1;
B3 = 0.1*[1 2 0 -2 -1];
T3 = tf(B3,A3);

% cascade (23T delay)
T = T1*T2*T3;
Aa = fliplr(T.den{1});
Ba = T.num{1};

% integration (18T delay)
Ws = fix(0.15*Fs);
B4 = ones(1,Ws)/Ws;
T4 = tf(B4,1);

% smoothing (4T delay)
B5 = ones(1,9)/9;
T5 = tf(B5,1);

% smoothing (1T delay)
%B5 = ones(1,3)/3;
%T5 = tf(B5,1);

% cascade (22T delay)
T = T4*T5;
Ab = fliplr(T.den{1});
Bb = T.num{1};

% apply filters
Signal = Signal - mean(Signal);
Signal = filter(Ba, Aa, Signal);
Result = filter(Bb, Ab, Signal.^2);
Delay = 23 + 22;