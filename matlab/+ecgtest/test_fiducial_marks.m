import ecgutilities.*

% le uma das derivaçoes do ecg
signal = interpret(x100,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay1] = mex.tompkins_filter(data, signal.fs);
R = mex.detect_qrs(sigI, signal.fs) - floor(delay1);
RR = diff(R);
RR = [RR(1); RR];

% detecta os pontos fiduciais
[sigD,sigM,gainD,delay2] = mex.yansun_filter(data, signal.fs);
R = R + floor(delay2);
F = fiducial_marks(sigD,sigM,R,RR,signal.fs);
data = [zeros(floor(delay2),1); data(1:end-floor(delay2))];

% visualizaçao
plot_fiducial_marks(data,F);
plot_fiducial_marks(sigD,F);
plot_fiducial_marks(sigM,F);