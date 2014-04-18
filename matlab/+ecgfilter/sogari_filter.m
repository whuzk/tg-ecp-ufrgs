function Result = sogari_filter(Signal,Fs,varargin)

if nargin > 2
    M = varargin{1};
else
    M = 3;
end

% linear filtering
Temp1 = smooth_and_diff(Signal);
%figure, plot(Temp1);

% non-linear transformation
Temp2 = ecgfilter.mobd(Temp1,M);
%figure, plot(Temp2);

% integration
Result = integrate(Temp2,Fs);


function Result = smooth_and_diff(Signal)
h = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
Result = wconv1(Signal, h, 'same');

function Result = integrate(Signal,Fs)
Ws = round(0.05*Fs)*2+1;
h = ones(1,Ws);
Result = wconv1(Signal,h,'same');