function Result = hermite(N, M, b)

% calculate vectors
Tb = math.tridiagonal_matrix(N, b);
[V,~] = eigs(Tb,M,'LA');

% correct polarity
Result = zeros(N,M);
k = floor(N/2)+1;
for m = 0:M-1
    Result(:,m+1) = V(:,m+1) .* sign(V(k,m+1)) * (-1)^floor(m/2);
end

function Result = tridiagonal_matrix(n, b)
i = 0:n-1;
j = 1:n-1;
tau = 1/(n*b^2);
D0 = -2*cos(pi*n*tau)*sin(pi*i*tau).*sin(pi*(n-i-1)*tau);
D1 = sin(pi*j*tau).*sin(pi*(n-j)*tau);
Result = diag(D0) + diag(D1,1) + diag(D1,-1);