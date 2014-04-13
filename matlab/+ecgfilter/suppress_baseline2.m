function Result = suppress_baseline2(Signal,Fs)

M = round(design(Fs,0.5,0.99));
len1 = 2*M+1;
len2 = 4*M+1;
B1 = 0.01*triangularPulse(1,len1,1:len1);
B2 = 0.02*triangularPulse(1,len2,1:len2);
Result = Signal - ecgmath.mmo_average(Signal,B1,B2);

function M = design(fs,f0,att)
M = fs*acos(att)/(2*pi*f0);