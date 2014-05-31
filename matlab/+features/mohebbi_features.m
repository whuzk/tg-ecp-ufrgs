function [C,STb] = mohebbi_features(signal, F, Beats, Template)

Fs = signal.fs;
gain = signal.gain;

% extraçao das caracteristicas
defJ = F(5,:)';
defTempJ = round(mean(F(5,1:30)));
[C,STb] = mex.mohebbi_features(Beats, Template, defJ, defTempJ, Fs, gain);