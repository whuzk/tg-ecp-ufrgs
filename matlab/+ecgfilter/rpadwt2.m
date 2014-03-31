function [A,D] = rpadwt2(x,J,H,G)

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
    for j = 1:J
        len = floor((length(A{j})+length(H)-1)/2);
        for p = 1:(len-length(A{j+1}))
            [A,D] = post(j,A,D,H,G);
        end
    end
elseif mod(i,2) == 0
    k = i/2;
    m = 1:min(i,length(H));
    A{j+1}(k) = A{j}(i-m+1)*H(m);
    D{j+1}(k) = A{j}(i-m+1)*G(m);
    [A,D] = rdwt(k,j+1,J,A,D,H,G);
end

function [A,D] = post(j,A,D,H,G)
k = length(A{j+1})+1;
i = 2*k;
d = max(0,i-length(A{j}));
m = 1:min(i,length(H)-d);
A{j+1}(k) = A{j}(i-d-m+1)*H(d+m);
D{j+1}(k) = A{j}(i-d-m+1)*G(d+m);