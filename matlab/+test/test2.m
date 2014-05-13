% le uma das derivaçoes do ecg
signal = utils.interpret(e0116,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay] = mex.tompkins_filter(data, signal.fs);
R = mex.detect_qrs(sigI, signal.fs);

% detecta os pontos fiduciais
[sigD,sigM] = mex.yansun_filter(data, signal.fs, delay);
F = mex.fiducial_marks(sigD, sigM, R, signal.fs);

% visualizaçao
data = [zeros(ceil(delay),1); data(1:end-ceil(delay))];
utils.plot_fiducial_marks(data,F);
utils.plot_fiducial_marks(sigD,F);
utils.plot_fiducial_marks(sigM,F);