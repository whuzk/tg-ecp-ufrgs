function [SignalH,SignalI,DelayH,DelayI,NewFs] = tompkins_preprocess(Signal, Fs)

% validate input
if ~isvector(Signal) || isempty(Signal)
  error('ecg must be a non-null vector');
end

% initialization
DelayH = 0;
DelayI = 0;
NewFs = 200;
resamp = (Fs ~= NewFs);

%% filter design

% resampling filter
if resamp
    [p,q] = rat(NewFs/Fs, 1E-12);
    pqmax = max(p,q);
    fc = 1/2/pqmax;
    L = 2*10*pqmax+1;
    f_r = [0 2*fc 2*fc 1];
    a_r = [1 1 0 0];
    h_r = firls(L-1, f_r, a_r);
    h_r = p*h_r.*kaiser(L,5)';
    d_r = floor(L/2/q);
end

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];
h_l = filter(b_l, a_l, [1 zeros(1,12)]);
d_l = 6;

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];
h_h = filter(b_h, a_h, [1 zeros(1,32)]);
d_h = 16;

% derivative filter
h_d = 1/8*[-1 -2 0 2 1];
d_d = 2;

% integration filter
Ws = round(0.15*NewFs);
h_i = ones(1,Ws)/Ws;
d_i = floor(Ws/2);

%% filtering

% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% resample signal if necessary
if resamp
    N = length(Signal);
    SignalUp(1:p:N*p) = Signal; % upsample
    Signal = conv2(SignalUp(:), h_r(:), 'full');
    Signal = Signal(1:q:end);   % downsample
    DelayH = DelayH + d_r;
end

% apply filters and update delays
SignalL = conv2(Signal, h_l(:), 'full');
SignalH = conv2(SignalL, h_h(:), 'full');
DelayH = DelayH + d_l + d_h;
SignalD = conv2(SignalH, h_d(:), 'full');
SignalI = conv2(SignalD.^2, h_i(:), 'full');
DelayI = DelayI + d_d + d_i;