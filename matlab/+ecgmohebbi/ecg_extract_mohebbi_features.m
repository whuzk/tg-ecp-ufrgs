function [F, D] = ecg_extract_mohebbi_features(Signal, Fs, R)
%   Extrai as caracteristicas do ECG pelo método de Mohebbi.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - vetor (K,1) indicando os picos de onda R
%
% Saída:
%   F - matriz (M,l) com l = 0.8*Fs carateristicas para cada uma das M
%       batidas cardiacas extraidas do ECG
%   D - matriz (K,1) indicando quais M das K batidas foram detectadas
%
import ecgmohebbi.*;
import utilities.*;

D = false(size(R));
N = length(Signal);
RefLen = fix(30*Fs);
Rlimit = fix(0.12*Fs);
Roffset = fix(0.12*Fs);
STlen = fix(0.16*Fs);

%% Parte I
R1 = R(R <= RefLen & R <= N-Rlimit-STlen+1);
if isempty(R1)
    F = [];
    return;
end

rows = ecg_extract_artifact_beats(Signal, Fs, R1);
R1(rows) = [];
if isempty(R1)
    F = [];
    return;
end

Template = ecg_construct_normal_template(Signal, Fs, R1 + Roffset);
Rtemplate = fix(length(Template)/2) - Roffset + 1;
Jtemplate = ecg_detect_jay_points(Template, Fs, Rtemplate);
TemplateST = Template(Jtemplate:min(length(Template),Jtemplate+STlen-1));

%% Parte II
R2 = R(R > RefLen);
if isempty(R2)
    F = [];
    return;
end
index1 = find(R > RefLen);

rows = ecg_extract_artifact_beats_based_on_template(...
    Signal, R2 + Roffset, Template);
R2(rows) = [];
index2 = find(~rows);

J2 = ecg_detect_jay_points(Signal, Fs, R2);
a = (J2 <= 0 | J2 > length(Signal)-STlen+1);
R2(a) = [];
J2(a) = [];
index3 = ~a;

F = ecg_prepare_st_segments(Signal, J2, TemplateST);
D(index1(index2(index3))) = true;

%% Plotagem de graficos
%{
figure;
hold on; grid on;
plot(1:length(Template), Template);
plot(Rtemplate, Template(Rtemplate), 'kx');
plot(Jtemplate, Template(Jtemplate), 'kx');
title('Template');
ecg_plot(TemplateST, 'Template ST segment');

figure;
hold on; grid on;
plot(1:N, Signal);
plot(R2, Signal(R2), 'kx');
plot(J2, Signal(J2), 'kx');
title('R peaks and J points');
%}