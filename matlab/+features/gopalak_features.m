function [C,S] = gopalak_features(signal, RR, Beats)

% extraçao das caracteristicas
[C,S] = mex.gopalak_features(Beats, RR, signal.fs);