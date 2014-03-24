function [SignalH,SignalI,Fs] = tompkins_preprocess(Signal, Fs)

% validate input
if ~isvector(Signal) || isempty(Signal)
  error('ecg must be a non-null vector');
end

%% filter design
resamp = false;%(Fs ~= 200);

% resampling filter
if resamp
    [p,q] = rat(200/Fs, 1E-12);
    pqmax = max(p,q);
    fc = 1/2/pqmax;
    L = 2*10*pqmax+1;
    f_r = [0 2*fc 2*fc 1];
    a_r = [1 1 0 0];
    h_r = firls(L-1, f_r, a_r);
    h_r = p*h_r.*kaiser(L,5)';
    d_r = floor(L/2/q);
    Fs = 200;
end

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];
h_l = filter(b_l, a_l, [1 zeros(1,12)]);
%d_l = 5;

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];
h_h = filter(b_h, a_h, [1 zeros(1,32)]);
%d_h = 15;

% alternative bandpass filter
%fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,13,20,30,1,20,Fs);
fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,18,25,30,1,20,Fs);
Hd1 = design(fbp,'equiripple');
h_b = Hd1.Numerator;
%d_b = 20;

% derivative filter (2T delay)
h_d = 1/8*[-1 -2 0 2 1];
%d_d = 2;

% integration filter (15T delay)
Ws = round(0.15*Fs);
h_i = ones(1,Ws)/Ws;
%d_i = floor(Ws/2);

%% filtering

% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% resample signal if necessary
if resamp
    N = length(Signal);
    SignalUp(1:p:N*p) = Signal; % upsample
    Signal = conv2(SignalUp(:), h_r(:), 'same');
    Signal = Signal(1:q:end);   % downsample
end

% apply filters and update delays
%SignalL = conv2(Signal, h_l(:), 'same');
%SignalH = conv2(SignalL, h_h(:), 'same');
SignalH = conv2(Signal, h_b(:), 'same');
SignalD = conv2(SignalH, h_d(:), 'same');
SignalI = conv2(SignalD.^2, h_i(:), 'same');