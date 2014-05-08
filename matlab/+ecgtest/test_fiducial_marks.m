import ecgutilities.*

% le uma das derivaçoes do ecg
signal = interpret(EDB.e0103,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay] = mex.tompkins_filter(data, signal.fs);
R = mex.detect_qrs(sigI, signal.fs);
RR = diff(R);
RR = [RR(1); RR];

% detecta os pontos fiduciais
[sigD,sigM] = mex.yansun_filter(data, signal.fs, delay);
F = mex.fiducial_marks(sigD,sigM,R,RR,signal.fs);
data = [zeros(ceil(delay),1); data(1:end-ceil(delay))];

% visualizaçao
plot_fiducial_marks(data,F);
plot_fiducial_marks(sigD,F);
plot_fiducial_marks(sigM,F);