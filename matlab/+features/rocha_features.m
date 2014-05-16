function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;
gain = signal.gain;

% extra�ao do ponto J de acordo com Pang
Jp = mex.pang_jpoints(RR, center, Fs);

% extra�ao dos pontos isoeletrico e J de acordo com Mohebbi
[Ib,Jb] = mex.mohebbi_ijpoints(Beats, F(3,:)', F(5,:)', Fs, gain);

% extra�ao dos desvios de segmento ST
C0 = mex.rocha_stdev(Beats, Ib, Jp, Jb);

% extra�ao dos segmentos de batida (demorado)
tic;
[S1,S2] = features.rocha_segments(Beats, Ib, Jb, F(end,:)');
toc;

% extra�ao dos coeficientes de hermite
[C1,C2] = mex.rocha_hermite(S1, S2);

% composi�ao das caracteristicas
C = [C0; C1; C2];