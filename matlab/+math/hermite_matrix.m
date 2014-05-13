function Result = hermite_matrix(N, M, l)
%   Obtem uma matriz de Hermite com base nas primeiras M funções de
%   Hermite com largura N e escala l.
%
Result = zeros(N,M);
t = (1:N)-ceil(N/2);
for m = 0:M-1
    Result(:,m+1) = math.hermite_function(m, t, l);
end
