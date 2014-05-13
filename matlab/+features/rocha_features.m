function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;

% extra�ao do ponto J de acordo com Pang
Jp = features.pang_jpoints(F(:,4), RR, Fs);

% extra�ao dos pontos isoeletrico e J de acordo com Rocha (demorado)
tic;
[Ir,Jr] = features.rocha_ijpoints(Beats, center, Fs);
toc;

% extra�ao dos segmentos de batida (demorado)
tic;
[S1,S2] = features.rocha_segments(Beats, Ir, Jr, F(:,end));
toc;

% extra�ao dos desvios de segmento ST
C0 = features.rocha_stdev(Beats, Ir, Jp, Jr);

% extra�ao dos coeficientes de hermite
[C1,C2] = features.rocha_hermite(S1, S2);

% composi�ao das caracteristicas
C = [C0; C1; C2];