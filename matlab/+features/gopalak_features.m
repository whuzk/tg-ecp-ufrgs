function [C,S] = gopalak_features(signal, RR, Beats)

% extraçao dos segmentos de batida
S = mex.gopalak_segments(Beats, RR, signal.fs);

% extraçao dos coefficientes de hermite
C = mex.gopalak_hermite(S);