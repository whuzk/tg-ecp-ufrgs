function [Result,delay] = sogari_filter(Signal,Fs,varargin)

delay = 0;

% process inputs
M = process_args(varargin,nargin-2);

% linear filtering
[Temp1,d] = smooth_and_diff(Signal);
delay = delay + d;

% improve dynamic range
[Temp2,G,d] = improve_range(Temp1,Fs,2^3,2^(floor(31/M))-1);
delay = delay + d;

% non-linear transformation
[Temp3,d] = ecgfilter.filter_mobd(Temp2,M);
delay = delay + d;

% integration
[Result,d] = integrate(Temp3,Fs);
delay = delay + d;

% plots
%{
figure, plot([Temp1 G]);
figure, plot(Temp2);
figure, plot(Temp3);
%}

function [Result,delay] = smooth_and_diff(Signal)
h = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
Result = filter(h,1,Signal);
delay = (length(h)-1)/2;

function [Result,G,delay] = improve_range(Signal,Fs,ming,maxg)
Signal = abs(Signal);
delay = floor(0.05*Fs);
G = ecgmath.running_max(Signal,2*Fs,-Inf);
Temp = Signal(1:end-delay)./max(ming,G(delay+1:end));
Result = [zeros(delay,1); min(maxg,maxg*Temp)];

function [Result,delay] = integrate(Signal,Fs)
Ws = floor(0.05*Fs)*2+1;
h = ones(1,Ws)/Ws;
Result = filter(h,1,Signal);
delay = (Ws-1)/2;

function M = process_args(args,na)
M = 3;  % default M-factor
if na > 0
    M = args{1};
end