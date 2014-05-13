function Result = baseline_filter(Signal,Fs)
import ecgmmo.*

M = round(0.0225*Fs+4);
B1 = ones(2*M+1,1);
B2 = ones(4*M+1,1);

oc = imerode(imdilate(imdilate(imerode(Signal,B1),B1),B2),B2);
co = imdilate(imerode(imerode(imdilate(Signal,B1),B1),B2),B2);
Result = Signal - (oc+co)/2;