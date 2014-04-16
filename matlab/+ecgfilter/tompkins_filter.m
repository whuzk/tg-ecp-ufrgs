function [SignalF,SignalI] = tompkins_filter(Signal, Fs)

%% filter design

% bandpass filter
fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,12,20,30,1,30,Fs);
Hd1 = design(fbp,'equiripple');
h_b = Hd1.Numerator;
%d_b = (length(h_f)-1)/2; ~[0.0817*Fs-0.0263]

% derivative filter
h_d = 1/8*[-1 -2 0 2 1];
%d_d = (length(h_d)-1)/2; =[2]

% integration filter
Ws = floor(0.15*Fs);
h_i = ones(1,Ws)/Ws;
%d_i = (Ws-1)/2; ~[0.075*Fs-0.5]

%% filtering
Signal = Signal(:) - Signal(1);
SignalF = conv2(Signal, h_b(:), 'same');
SignalD = conv2(SignalF, h_d(:), 'same');
SignalI = conv2(SignalD.^2, h_i(:), 'same');