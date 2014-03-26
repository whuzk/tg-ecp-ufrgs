function [SignalF,SignalI,Fs] = tompkins_preprocess_original(Signal, Fs)

% validate input
if ~isvector(Signal) || isempty(Signal)
  error('ecg must be a non-null vector');
end

%% filter design
resamp = (Fs ~= 200);

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
    %d_r = floor(L/2/q);
    Fs = 200;
end

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];

% derivative filter
h_d = 1/8*[-1 -2 0 2 1];

% cascade of previous filters
t_f = minreal(tf(b_l,a_l)*tf(b_h,a_h)*tf(h_d,1));
a_f = fliplr(t_f.den{1});
b_f = t_f.num{1};
h_f = filter(b_f, a_f, [1; zeros(length(b_f)-1,1)]);
%d_f = floor(length(b_f)/2);

% integration filter
Ws = round(0.15*Fs);
h_i = ones(1,Ws)/Ws;
%d_i = floor(Ws/2);

%% filtering
% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% resample if necessary
if resamp
    N = length(Signal);
    SignalUp(1:p:N*p) = Signal; % upsample
    Signal = conv2(SignalUp(:), h_r(:), 'same');
    Signal = Signal(1:q:end);   % downsample
end

% apply filters
SignalF = conv2(Signal, h_f(:), 'same');
SignalI = conv2(SignalF.^2, h_i(:), 'same');