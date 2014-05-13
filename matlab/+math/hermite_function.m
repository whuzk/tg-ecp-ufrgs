function Result = hermite_function(n, t, l)
% Obtem os valores da funçao de hermite de ordem n, com escala l,
% avaliados nos pontos de entrada t.
x = t./l;
P = math.hermite_poly(n);
C = 1/sqrt(prod(1:n)*(2^n)*sqrt(pi)*l);
Result = C*exp(-0.5*x.^2).*polyval(P,x);