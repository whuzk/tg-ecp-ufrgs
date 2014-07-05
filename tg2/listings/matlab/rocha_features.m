function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;
gain = signal.gain;

% extracao do ponto J de acordo com Pang
Jp = features.pang_jpoints(F(4,:)', RR, Fs);

% extracao dos pontos isoeletrico e J de acordo com Rocha
[Ir,Jr] = features.rocha_ijpoints(Beats, center, F(3,:)', F(5,:)', Fs, gain);

% extracao dos segmentos de batida
[S1,S2] = features.rocha_segments(Beats, Ir, Jr, F(end,:)');

% extracao dos desvios de segmento ST
C0 = features.rocha_stdev(Beats, Ir, Jp, Jr);

% extracao dos coeficientes de hermite
[C1,C2] = features.rocha_hermite(S1, S2);

% composicao das caracteristicas
C = [C0; C1; C2];