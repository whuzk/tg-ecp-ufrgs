function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;
gain = signal.gain;

% extraçao do ponto J de acordo com Pang
Jp = mex.pang_jpoints(RR, center, Fs);

% extraçao dos pontos isoeletrico e J de acordo com Mohebbi
[Ib,Jb] = mex.mohebbi_ijpoints(Beats, F(3,:)', F(5,:)', Fs, gain);

% extraçao dos desvios de segmento ST
C0 = mex.rocha_stdev(Beats, Ib, Jp, Jb);

% extraçao dos segmentos de batida (demorado)
tic;
[S1,S2] = features.rocha_segments(Beats, Ib, Jb, F(end,:)');
toc;

% extraçao dos coeficientes de hermite
[C1,C2] = mex.rocha_hermite(S1, S2);

% composiçao das caracteristicas
C = [C0; C1; C2];