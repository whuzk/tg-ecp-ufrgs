function [R,RR,F,Beats,Templates,delay] = preprocess(signal)

% obtem as amostras do sinal e a taxa de amostragem
data = signal.data - signal.inival;
Fs = signal.fs;

% detecta os picos de onda R
[sigI,delay] = mex.tompkins_filter(data, Fs);
[R,RR] = mex.detect_qrs_prod(sigI, Fs);

% detecta os pontos fiduciais
[sigD,sigM] = mex.yansun_filter(data, Fs, delay);
F = mex.fiducial_marks(sigD, sigM, R, Fs);

% extrai as batidas
sigN = mex.noise_filter(data, Fs, delay);
[Beats,F] = mex.extract_beats(sigN, F, Fs);

% constroi versoes de template e detecta artefatos
[Templates,idx] = mex.template_and_artifacts(Beats, 30);

% remove os artefatos
R(idx) = [];
RR(idx) = [];
F(idx,:) = [];
Beats(:,idx) = [];
Templates(:,idx) = [];