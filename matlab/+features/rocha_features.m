function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;

% extraçao do ponto J de acordo com Pang
Jp = features.pang_jpoints(F(:,4), RR, Fs);

% extraçao dos pontos isoeletrico e J de acordo com Rocha (demorado)
tic;
[Ir,Jr] = features.rocha_ijpoints(Beats, center, Fs);
toc;

% extraçao dos segmentos de batida (demorado)
tic;
[S1,S2] = features.rocha_segments(Beats, Ir, Jr, F(:,end));
toc;

% extraçao dos desvios de segmento ST
C0 = features.rocha_stdev(Beats, Ir, Jp, Jr);

% extraçao dos coeficientes de hermite
[C1,C2] = features.rocha_hermite(S1, S2);

% composiçao das caracteristicas
C = [C0; C1; C2];