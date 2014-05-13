function [C,STb] = mohebbi_features(signal, F, Beats, Template)

center = floor(size(Beats,1)/2)+1;
defJ = mean(F(1:30,5));
Fs = signal.fs;
gain = signal.gain;

% extra�ao dos pontos isoeletrico e J de acordo com Mohebbi
Jt = features.mohebbi_jpoints(Template, center, defJ, Fs, gain);
Jb = features.mohebbi_jpoints(Beats, center, F(:,5), Fs, gain);

% extra�ao dos segmentos ST (demorado)
STt = features.mohebbi_segments(Template, Jt, Fs);
STb = features.mohebbi_segments(Beats, Jb, Fs);

% extra�ao das diferen�as entre os segmentos e o do template
C = features.mohebbi_stdiff(STt, STb);