function [C,S] = gopalak_features(RR, Beats)

% extraçao dos segmentos de batida (demorado)
S = features.gopalak_segments(Beats, RR);

% extraçao dos coefficientes de hermite
C = features.gopalak_hermite(S);