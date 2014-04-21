function Result = sogari_conv(Signal,Fs,varargin)

% process inputs
[M,B] = process_args(varargin,nargin-2);

% linear filtering
Temp1 = smooth_and_diff(Signal);

% absolute value
Temp2 = abs(Temp1);

% improve dynamic range
[Temp3,G] = improve_range(Temp2,Fs,B);

% non-linear transformation
Temp4 = ecgfilter.conv_mobd(Temp3,M);

% integration
Result = integrate(Temp4,Fs);

% plots
%{
figure, plot(Temp1);
figure, plot([Temp2 G]);
figure, plot(Temp3);
figure, plot(Temp4);
%}

function Result = smooth_and_diff(Signal)
h = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
Result = wconv1(Signal, h, 'same');

function [Result,G] = improve_range(Signal,Fs,B)
maxq = 2^B-1;
G = ecgmath.running_max(Signal,2*Fs,maxq);
Result = min(maxq,maxq*Signal./G);

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