% le uma das derivaçoes do ecg
recordIndex = 2;
signal = utils.interpret_ecg(x100, recordIndex);
data = signal.data - signal.inival;
Fs = signal.fs;

% filtragem e detecçao
[sigI,sigD,sigM,sigN,delay] = mex.preprocess_filter(data, Fs);
R = mex.detect_qrs(sigI, Fs);
F = mex.fiducial_marks(sigD, sigM, R, Fs);

% visualizaçao
data = [zeros(ceil(delay),1); data(1:end-ceil(delay))];
plot.plot_fiducial_marks(data,F);
plot.plot_fiducial_marks(sigD,F);
plot.plot_fiducial_marks(sigM,F);