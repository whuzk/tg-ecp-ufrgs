%% ECG data
signal = utils.interpret_ecg(e0103, 1);
data = int16(signal.data - signal.inival);
simin = timeseries(data, signal.time);
Fs = signal.fs;
Fm = 50;

%% QRS filters
[qrsbl,qrsal,qrsgl,qrsdl] = intfdesign.lowpass('N,F3db',2,11,Fs);
[qrsbh,qrsah,qrsgh,qrsdh] = intfdesign.highpass('N,F3db',1,5,Fs);
[qrsbd,qrsad,qrsgd,qrsdd] = intfdesign.derivative('N,M',1,3);
[qrsbi,qrsai,qrsgi,qrsdi] = intfdesign.maverage('Width',0.15,Fs);
delay_qrs = qrsdl + qrsdh + qrsdd + qrsdi;

%% FP filters
[fpbl,fpal,fpgl,fpdl] = intfdesign.lowpass('N,F3db',2,11,Fs);
[fpbi,fpai,fpgi,fpdi] = intfdesign.maverage('Width',0.05,Fs);
[fpbd,fpad,fpgd,fpdd] = intfdesign.derivative('N,M',1,0);
s = round(0.06*Fs); fpeds = 2*s+1;
delay_fp1 = ceil(delay_qrs-fpdl-fpdi-fpdd);
delay_fp2 = ceil(delay_qrs-fpdl-s);

%% Noise filters
[noisebs,noiseas] = cheby1(2,0.1,[Fm-1 Fm+1]*2/Fs,'stop');
noiseds = mean(grpdelay(noisebs,noiseas,[5 15]*2/Fs));
[noisebl,noiseal] = butter(4,40*2/Fs);
noisedl = mean(grpdelay(noisebl,noiseal,[5 15]*2/Fs));
delay_noise = ceil(delay_qrs-noiseds-noisedl);

%% Detection
BufLen = 4*Fs;
L1 = floor(0.10*Fs);
L2 = floor(0.02*Fs);
L3 = floor(0.15*Fs);
L4 = floor(0.08*Fs);
FrameLen = 2*floor(0.6*Fs)+1;
TempCount = 30;

%% Extraction
mawidth = floor(0.05*Fs);
beatlpb = [1 zeros(1,mawidth-1) -1];
beatlpa = [1 -1];
beatlpd = (mawidth-1)/2;
beatdeb = [1 -1];
beatdea = 1;
beatded = 0.5;
stsegsize = floor(0.08*Fs)*2;
jaythresh = floor(mawidth*signal.gain*1.25/Fs);
newL1 = floor(0.02*Fs);
newL2 = floor(0.12*Fs);