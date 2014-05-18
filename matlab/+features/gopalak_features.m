function [C,S] = gopalak_features(signal, RR, Beats)

% extra�ao dos segmentos de batida
S = mex.gopalak_segments(Beats, RR, signal.fs);

% extra�ao dos coefficientes de hermite
C = mex.gopalak_hermite(S);