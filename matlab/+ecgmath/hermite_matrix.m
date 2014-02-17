function Result = hermite_matrix(w, m, l)
%   Obtem uma matriz de Hermite com base nas primeiras m funções de
%   Hermite de largura w e escala l.
%
% Entradas:
%   w - largura do dominio das funçoes de hermite
%   m - numero de funçoes a serem consideradas
%   l - fator de escala das funçoes de hermite
%
% Saída:
%   matriz (m,w) de regressao hermitiana
%
Result = zeros(w,m);
X = (1:w)-floor(w/2);
for i = 1:m
    Result(:,i) = ecgmath.hermite_function(X, i-1, l);
end
