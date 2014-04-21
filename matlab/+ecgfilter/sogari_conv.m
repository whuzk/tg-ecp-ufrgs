function Result = sogari_conv(Signal,Fs,varargin)

% process inputs
[M,B] = process_args(varargin,nargin-2);

% linear filtering
Temp1 = smooth_and_diff(Signal);

% improve dynamic range
[Temp2,G] = improve_range(Temp1,Fs,2^3,2^(floor(31/M))-1);

% non-linear transformation
Temp3 = ecgfilter.conv_mobd(Temp2,M);

% integration
Result = integrate(Temp3,Fs);

% plots
%{
figure, plot([Temp1 G]);
figure, plot(Temp2);
figure, plot(Temp3);
%}

function Result = smooth_and_diff(Signal)
h = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
Result = wconv1(Signal, h, 'same');

function [Result,G] = improve_range(Signal,Fs,ming,maxg)
Signal = abs(Signal);
G = ecgmath.running_max(Signal,2*Fs,-Inf);
Result = min(maxg,maxg*Signal./max(ming,G));

function Result = integrate(Signal,Fs)
Ws = floor(0.05*Fs)*2+1;
h = ones(1,Ws)/Ws;
Result = wconv1(Signal,h,'same');

function [M,B] = process_args(args,na)
M = 3;  % default M-factor
B = 8;  % default resolution
if na > 0
    M = args{1};
end
if na > 1
    B = args{2};
end