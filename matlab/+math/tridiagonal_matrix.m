function Result = tridiagonal_matrix(n, b)
% Obtem a matriz tridiagonal quadrada que comuta com a matriz de Fourier
% centralizada e possui dilataçao especificada pelo parametro b.
%
i = 0:n-1;
j = 1:n-1;
tau = 1/(n*b^2);
D0 = -2*cos(pi*n*tau)*sin(pi*i*tau).*sin(pi*(n-i-1)*tau);
D1 = sin(pi*j*tau).*sin(pi*(n-j)*tau);
Result = diag(D0) + diag(D1,1) + diag(D1,-1);