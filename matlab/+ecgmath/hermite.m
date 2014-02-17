function Result = hermite(n, l, b)
%   Obtem as primeiras l fun�oes de Hermite discretas de tamanho n com
%   parametro de dilata�ao b. Calculo feito atraves dos autovetores e
%   autovalores de uma matriz tridiagonal, de acordo com o metodo proposto
%   por Mugler et al.
%
% Entradas:
%   n - tamanho da matriz
%   l - numero de fun�oes de Hermite
%   b - parametro de dilata�ao
%
% Sa�da:
%   matriz (n,l) com as fun�oes de Hermite discretas
%
if (l > n)
    error('number of functions cannot exceed function size');
end

Tb = tridiagonal_matrix(n, b);
[Result,~] = eigs(Tb,l,'LA');

i = fix(n/2);
for j = 1:l
    Result(:,j) = Result(:,j) .* sign(Result(i,j));
end


function Result = tridiagonal_matrix(n, b)
%   Obtem a matriz tridiagonal quadrada que comuta com a matriz de Fourier
%   centralizada e possui dilata�ao especificada pelo parametro b.
%
% Entradas:
%   n - tamanho da matriz
%   b - parametro de dilata�ao
%
% Sa�da:
%   matriz tridiagonal
%
i = 0:n-1;
j = 1:n-1;
tau = 1/(n*b^2);
D0 = -2*cos(pi*n*tau)*sin(pi*i*tau).*sin(pi*(n-i-1)*tau);
D1 = sin(pi*j*tau).*sin(pi*(n-j)*tau);
Result = diag(D0) + diag(D1,1) + diag(D1,-1);
