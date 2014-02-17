function Result = hermite_function(X, n, l)
%   Obtem os valores da funçao de hermite de ordem n, com escala l,
%   avaliados nos pontos de entrada X.
%
% Entradas:
%   X - vetor de entrada da funçao
%   n - ordem da funçao de hermite
%   l - fator de escala da funçao
%
% Saída:
%   valores da funçao de hermite avaliados nos pontos de entrada
%
P = hermite_poly(n);
C = 1/sqrt(factorial(n)*(2^n)*sqrt(pi)*l);
Result = C*exp(-0.5*(X/l).^2).*polyval(P,X);

function Result = hermite_poly(n)
%   Obtem o polinomio de hermite de ordem n.
%
% Entradas:
%   n - ordem do polinomio de hermite
%
% Saída:
%   vetor (1,n+1) com os coeficientes do polinomio de hermite
%
Result = zeros(1,n+1);
Result(1) = 2^n;

for i = 1:floor(n/2)
    if 2*i == n
        CUM = cumprod(i+1:n);
        Result(end) = (-1)^i*CUM(end);
    else
        NUM = cumprod(i+1:n);
        DEN = cumprod(1:n-2*i);
        Result(2*i+1) = (-1)^i * 2^(n-2*i) * NUM(end) / DEN(end);
    end
end
