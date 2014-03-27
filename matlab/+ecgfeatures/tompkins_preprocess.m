function [SignalF,SignalI,Fs] = tompkins_preprocess(Signal, Fs)

% validate input
if ~isvector(Signal) || isempty(Signal)
  error('ecg must be a non-null vector');
end

%% filter design

% bandpass filter
%fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,13,20,30,1,20,Fs);
fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,18,25,30,1,20,Fs);
Hd1 = design(fbp,'equiripple');
h_f = Hd1.Numerator;
%d_f = floor(length(h_f)/2);

% derivative filter
h_d = 1/8*[-1 -2 0 2 1];
%d_d = floor(length(h_d)/2);

% integration filter
Ws = round(0.15*Fs);
h_i = ones(1,Ws)/Ws;
%d_i = floor(Ws/2);

%% filtering
% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% apply filters
SignalF = conv2(Signal, h_f(:), 'same');
SignalD = conv2(SignalF, h_d(:), 'same');
SignalI = conv2(SignalD.^2, h_i(:), 'same');