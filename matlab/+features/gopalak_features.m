function [C,S] = gopalak_features(RR, Beats)

% extra�ao dos segmentos de batida (demorado)
S = features.gopalak_segments(Beats, RR);

% extra�ao dos coefficientes de hermite
C = features.gopalak_hermite(S);