function [A,D] = rpadwt(x,J,H,G)

N = length(x);
A = cell(J+1,1);
D = cell(J+1,1);
A{1} = x(:)';
H = H(:);
G = G(:);
for i = 1:N
    [A,D] = rdwt(i,1,J,A,D,H,G);
end
A = A(2:end);
D = D(2:end);

function [A,D] = rdwt(i,j,J,A,D,H,G)
if j > J
    return;
elseif mod(i,2) ~= 0
    k = (i+1)/2;
    L = min(i,length(H));
    m = 1:L;
    A{j+1}(k) = A{j}(i-m+1)*H(m);
    D{j+1}(k) = A{j}(i-m+1)*G(m);
else
    [A,D] = rdwt(i/2,j+1,J,A,D,H,G);
end