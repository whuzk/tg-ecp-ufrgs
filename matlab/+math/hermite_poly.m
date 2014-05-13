function Result = hermite_poly(n)
% Obtem o polinomio de hermite de ordem n.
Result = zeros(1,n+1);
for i = 0:floor(n/2)
    k = 2*i;
    Result(k+1) = (-1)^i * 2^(n-k) * prod(i+1:n) / prod(1:n-k);
end