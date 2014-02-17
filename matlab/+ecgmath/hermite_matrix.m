function Result = hermite_matrix(w, m, l)
%   Obtem uma matriz de Hermite com base nas primeiras m fun��es de
%   Hermite de largura w e escala l.
%
% Entradas:
%   w - largura do dominio das fun�oes de hermite
%   m - numero de fun�oes a serem consideradas
%   l - fator de escala das fun�oes de hermite
%
% Sa�da:
%   matriz (m,w) de regressao hermitiana
%
Result = zeros(w,m);
X = (1:w)-floor(w/2);
for i = 1:m
    Result(:,i) = ecgmath.hermite_function(X, i-1, l);
end
