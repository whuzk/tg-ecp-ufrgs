function [C,S] = gopalak_features(RR, Beats)

% extracao dos segmentos de batida
S = features.gopalak_segments(Beats, RR);

% extracao dos coefficientes de hermite
C = features.gopalak_hermite(S);