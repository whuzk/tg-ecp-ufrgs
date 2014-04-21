function Result = trahanias_filter(Signal,Fs)
import ecgmmo.*

% noise suppression
B = [0 1 5 1 0];
oc = mmo_open_closing(Signal,B,B);
co = mmo_close_opening(Signal,B,B);
Temp1 = (oc+co)/2;

% peak-valley-extractor
B1 = zeros(1,13);
Temp2 = mmo_open_closing(Temp1,B1,B1);

% result signal
Result = Temp1 - Temp2;