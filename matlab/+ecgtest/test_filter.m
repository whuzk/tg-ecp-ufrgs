[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);

% low-pass filter
b_l = 1/36*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_l = [1 -2 1];
%d_l = 5;

% high-pass filter
b_h = 1/32*[-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_h = [1 -1];
%d_h = 15;

% alternative bandpass filter
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

% cascade of filters
t1 = minreal(tf(b_l,a_l)*tf(b_h,a_h)*tf(h_d,1));
a1 = fliplr(t1.den{1});
b1 = t1.num{1};
%d1 = floor(length(b1)/2);
h1 = filter(b1, a1, [1; zeros(length(b1)-1,1)]);
h2 = conv(h_b, h_d);

% apply filters
SignalF1 = conv2(Signal, h1(:), 'same').^2;
SignalI1 = conv2(SignalF1, h_i(:), 'same');

SignalF2 = conv2(Signal, h2(:), 'same').^2;
SignalI2 = conv2(SignalF2, h_i(:), 'same');

figure, plot([SignalI2 SignalI1]);