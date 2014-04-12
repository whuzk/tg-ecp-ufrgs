function Result = suppress_noise(Signal, Fs)

flp = fdesign.lowpass('Fp,Fst,Ap,Ast',30,50,0.01,60,Fs);
hd = design(flp,'equiripple');
Result = conv2(Signal, hd.numerator(:), 'same');