function Result = ecg_classify_ischemic_beats(F)
%   Classifica as batidas de acordo com a eleva��o e depress�o do segmento
%   ST, atraves de uma rede neural tendo como entrada as caracteristicas
%   extraidas de cada batida. 
%
% Entradas:
%   F - conjunto de caracteristicas das batidas
%
% Sa�da:
%   indice das batidas classificadas como isquemicas
%
global RochaNet;

O = sim(RochaNet, F')';
Result = O(:,1) > 0 | O(:,2) > 0;
