import ecgutilities.*

% le uma das derivaçoes do ecg
signal = interpret(e0103,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay1] = mex.tompkins_filter(data, signal.fs);
R = mex.detect_qrs(sigI, signal.fs) - floor(delay1);
R = adjust_qrs(data,signal.lead,R,signal.fs);

% detecta os pontos fiduciais
[sigD,sigM,gainD,delay2] = mex.yansun_filter(data, signal.fs);
R = R + floor(delay2);
RR = diff([0; R]);
F = utils.fiducial_marks(sigD,gainD,sigM,signal.lead,R,RR,signal.fs);
data = [zeros(floor(delay2),1); data(1:end-floor(delay2))];

% visualizaçao
plot_fiducial_marks(data,R,F);
plot_fiducial_marks(sigD,R,F);
plot_fiducial_marks(sigM,R,F);