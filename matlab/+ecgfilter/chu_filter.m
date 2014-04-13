function Result = chu_filter(Signal,Fs)
import ecgmath.*

% noise suppression
B = 0.01*[0 1 5 1 0];
oc = mmo_open_closing(Signal,B,B);
co = mmo_close_opening(Signal,B,B);
Temp1 = (oc+co)/2;

% baseline removal
M = round(design(Fs,0.5,0.9));
len1 = 2*M+1;
len2 = 4*M+1;
B1 = 0.01*triangularPulse(1,len1,1:len1);
B2 = 0.02*triangularPulse(1,len2,1:len2);
oc = mmo_open_closing(Temp1,B1,B2);
co = mmo_close_opening(Temp1,B1,B2);
Temp2 = (oc+co)/2;

% result signal
Result = Temp1 - Temp2;


function M = design(fs,f0,att)
M = fs*acos(att)/(2*pi*f0);