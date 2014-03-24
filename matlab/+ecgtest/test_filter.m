[Fs,~,~,~,Data] = ecgutilities.interpret(Database2.e0108);
Signal = Data(:,2);

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

% apply filters
SignalL1 = conv2(Signal, h_l(:), 'same');
SignalH1 = conv2(SignalL1, h_h(:), 'same');
SignalD1 = conv2(SignalH1, h_d(:), 'same');
SignalI1 = conv2(SignalD1.^2, h_i(:), 'same');

SignalH2 = conv2(Signal, h_b(:), 'same');
SignalD2 = conv2(SignalH2, h_d(:), 'same');
SignalI2 = conv2(SignalD2.^2, h_i(:), 'same');

figure, plot([SignalH2 SignalH1]);
figure, plot([SignalI2 SignalI1]);