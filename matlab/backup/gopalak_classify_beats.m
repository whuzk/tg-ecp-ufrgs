function Result = ecg_classify_ischemic_beats(F)
%   Classifica as batidas em isquemicas ou nao, atraves de uma rede neural
%   tendo como entrada as caracteristicas extraidas de cada batida.
%
% Entradas:
%   F - caracteristicas das batidas
%
% Saída:
%   indice das batidas classificadas como isquemicas
%
global GopalakNets;

k = length(GopalakNets);
m = size(F,1);
decision = zeros(m,1);
for i = 1:k
    O = sim(GopalakNets{i}, F')';
    C = O(:,1) > 0 | O(:,2) > 0;
    decision = decision + double(C);
end
Result = decision > k/2;
