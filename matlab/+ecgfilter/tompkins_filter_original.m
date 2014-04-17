function [SignalF,SignalI] = tompkins_filter_original(Signal, Fs)

%% filter design

% lowpass
b_l = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];

% highpass
b_h = [-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32];
a_h = [1 -1];

% cascade
t_f = minreal(tf(b_l,a_l)*tf(b_h,a_h));
a_f = fliplr(t_f.den{1});
b_f = t_f.num{1};
h_f = filter(b_f, a_f, [1 zeros(1,length(b_f))]);
%d_f = (length(h_f)-1)/2;

% derivative
h_d = 1/8*[2 1 0 -1 -2];
%d_d = (length(h_d)-1)/2; =[2]

% integration
h_i = ones(1,32)/32;
%d_i = (length(h_i)-1)/2;

%% filtering
Signal = Signal(:) - Signal(1);
SignalF = wconv1(Signal, h_f, 'same');
SignalD = wconv1(SignalF, h_d, 'same');
SignalI = wconv1(SignalD.^2, h_i, 'same');