function Result = hermite(N, M, b)
%   Obtem as primeiras M funçoes de Hermite discretas de tamanho N com
%   parametro de dilataçao b. Calculo feito atraves dos autovetores e
%   autovalores de uma matriz tridiagonal, de acordo com o metodo proposto
%   por Mugler et al.
%
Tb = math.tridiagonal_matrix(N, b);
[V,~] = eigs(Tb,M,'LA');

% correct function polarity
Result = zeros(N,M);
k = floor(N/2)+1;
for m = 0:M-1
    Result(:,m+1) = V(:,m+1) .* sign(V(k,m+1)) * (-1)^floor(m/2);
end