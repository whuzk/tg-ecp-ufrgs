function [C,S] = gopalak_features(signal, RR, Beats)

% extra�ao das caracteristicas
[C,S] = mex.gopalak_features(Beats, RR, signal.fs);