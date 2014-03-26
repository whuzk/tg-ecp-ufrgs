[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];

% cascade of previous filters
t_f = minreal(tf(b_l,a_l)*tf(b_h,a_h));
a_f = fliplr(t_f.den{1});
b_f = t_f.num{1};
h_f = filter(b_f, a_f, [1; zeros(length(b_f)-1,1)]);

% alternative bandpass filter
fbp = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1,8,18,25,30,1,20,Fs);
Hd1 = design(fbp,'equiripple');
h_b = Hd1.Numerator;

% derivative filter (2T delay)
h_d = 1/8*[-1 -2 0 2 1];

% integration filter (15T delay)
Ws = round(0.15*Fs);
h_i = ones(1,Ws)/Ws;

% apply filters
SignalF1 = conv2(Signal, h_f(:), 'same');
SignalD1 = conv2(SignalF1, h_d(:), 'same');
SignalI1 = conv2(SignalD1.^2, h_i(:), 'same');

SignalF2 = conv2(Signal, h_b(:), 'same');
SignalD2 = conv2(SignalF2, h_d(:), 'same');
SignalI2 = conv2(SignalD2.^2, h_i(:), 'same');

% plot figures
figure, plot([SignalI2 SignalI1]);