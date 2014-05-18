function [C,STb] = mohebbi_features(signal, F, Beats, Template)

defI = round(mean(F(3,1:30)));
defJ = round(mean(F(5,1:30)));
Fs = signal.fs;
gain = signal.gain;

% extra�ao dos pontos isoeletrico e J de acordo com Mohebbi
[~,Jt] = mex.mohebbi_ijpoints(Template, defI, defJ, Fs, gain);
[~,Jb] = mex.mohebbi_ijpoints(Beats, F(3,:)', F(5,:)', Fs, gain);

% extra�ao dos segmentos ST
STt = mex.mohebbi_segments(Template, Jt, Fs);
STb = mex.mohebbi_segments(Beats, Jb, Fs);

% extra�ao das diferen�as entre os segmentos e o do template
C = mex.mohebbi_stdiff(STt, STb);