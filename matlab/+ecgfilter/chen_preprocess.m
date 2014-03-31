function Result = chen_preprocess(Signal, Fs, W)

% validate input
if ~isvector(Signal) || isempty(Signal)
  error('ecg must be a non-null vector');
end

% vectorize and eliminate drift
Signal = Signal(:) - Signal(1);

% wavelet de-noise
SignalW = ecgfilter.rtdenoise(Signal,3,256,W);
%[b,a] = butter(4,30*2/Fs);
%SignalW = filter(b,a,Signal);

% linear highpass filter
SignalH = lhpf(SignalW,5);

% non-linear lowpass filter
Result = nlpf(SignalH,Fs);

% plot
figure, plot(SignalW);
figure, plot(SignalH);


function Result = lhpf(Signal,M)
h1 = ones(1,M)/M;
h2 = [zeros(1,(M-1)/2) 1 zeros(1,(M-1)/2)];
h = h2-h1;
Result = wconv1(Signal,h,'same');

function Result = nlpf(Signal,Fs)
Wn = round(0.15*Fs);
h = ones(1,Wn)/Wn;
Result = wconv1(Signal.^2,h,'same');