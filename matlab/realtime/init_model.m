%% ECG data
signal = utils.interpret_ecg(e0103, 1);
data = signal.data - signal.inival;
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
delay_noise = ceil(noiseds + noisedl);