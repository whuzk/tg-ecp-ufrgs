function Result = ecg_classify_ischemic_beats(F)
%   Classifica as batidas em isquemicas ou nao, atraves de uma rede neural
%   tendo como entrada as caracteristicas extraidas de cada batida.
%
% Entradas:
%   F - dados que caracterizam o segmento ST das batidas
%
% Saída:
%   indice das batidas classificadas como isquemicas
%
global MohebbiNet;

O = sim(MohebbiNet, F')';
Result = O(:,1) > 0 & O(:,2) < 0;
