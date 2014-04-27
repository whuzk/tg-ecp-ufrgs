function [Result,delay] = suppress_noise(Signal, Fs)

flp = fdesign.lowpass('Fp,Fst,Ap,Ast',30,50,0.01,60,Fs);
hd = design(flp,'equiripple');
h = hd.numerator(:);
Result = filter(h, 1, Signal);
delay = (length(h)-1)/2;