function Result = krishnan_filter(Signal,Fs)
import ecgmmo.*

% baseline removal
L1 = round(0.1*Fs)*2+1;
L2 = round(0.75*L1)*2+1;
B1 = zeros(1,L1);
B2 = zeros(1,L2);
Temp1 = mmo_open_closing(Signal,B1,B2);
Temp2 = Signal - Temp1;

% noise suppression
B1 = [0 1 5 1 0];
B2 = [0 0 0 0 0];
op = mmo_opening(Temp2,B1,B2);
cl = mmo_closing(Temp2,B1,B2);
Result = (op+cl)/2;