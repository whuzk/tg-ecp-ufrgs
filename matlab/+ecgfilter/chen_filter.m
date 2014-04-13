function Result = chen_filter(Signal, Fs, wname)

% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% wavelet de-noise
SignalW = wvlf(Signal,3,wname);

% linear highpass filter
SignalH = lhpf(SignalW,5);

% non-linear lowpass filter
Result = nlpf(SignalH,Fs);

% plot
%figure, plot(SignalW);
%figure, plot(SignalH);


function Y = wvlf(X,J,wname)
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);
[A,D] = ecgfilter.dpadwt2(X,J,Lo_D,Hi_D);
D = ecgfilter.wavelet_denoise(D,J,4);
Y = ecgfilter.dpaidwt2(A{end},D,length(X),Lo_R,Hi_R);

function Result = lhpf(Signal,M)
h1 = ones(1,M)/M;
h2 = [zeros(1,(M-1)/2) 1 zeros(1,(M-1)/2)];
Result = wconv1(Signal,h2-h1,'same');

function Result = nlpf(Signal,Fs)
Wn = round(0.15*Fs);
h = ones(1,Wn)/Wn;
Result = wconv1(Signal.^2,h,'same');