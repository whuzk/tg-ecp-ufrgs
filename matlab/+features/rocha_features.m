function [C,S1,S2] = rocha_features(signal, RR, F, Beats)

%center = floor(size(Beats,1)/2)+1;
Fs = signal.fs;
gain = signal.gain;
defI = F(3,:)';
defJ = F(5,:)';
endPts = F(end,:)';

% extraçao das caracteristicas
[C0,C1,C2,S1,S2] = ...
    mex.c_rocha_features(Beats, RR, defI, defJ, endPts, Fs, gain);

% composiçao das caracteristicas
C = [C0; C1; C2];