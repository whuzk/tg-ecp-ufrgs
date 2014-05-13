% le uma das derivaçoes do ecg
signal = utils.interpret_ecg(e0103,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay] = mex.tompkins_filter(data, signal.fs);
R = mex.detect_qrs(sigI, signal.fs);

% detecta os pontos fiduciais
[sigD,sigM] = mex.yansun_filter(data, signal.fs, delay);
F = mex.fiducial_marks(sigD, sigM, R, signal.fs);

% extrai as batidas
sigN = mex.noise_filter(data, signal.fs, delay);
[Beats,F] = mex.extract_beats(sigN, F, signal.fs);

% construcao de template e deteccao de artefatos
[Templates,idx] = mex.template_and_artifacts(Beats,30);

% visualizacao
figure, plot(Templates(:,30));