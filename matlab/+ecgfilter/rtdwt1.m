function [A,D] = rtdwt1(X,J,H,G)

N = length(X);
L = length(H);
A = cell(J+1,1);
D = cell(J+1,1);
S = zeros(J+1,L);
H = H(:);
G = G(:);

len = N;
for j = 1:J
    len = floor(len/2);
    A{j+1} = zeros(1,len);
    D{j+1} = zeros(1,len);
end

k = 1;
for i = 1:N
    S(1,:) = push(S(1,:), X(i));
    [S,A,D] = recursive(k,1,J,S,A,D,H,G);
    k = mod(k,2^J)+1;
end

A{1}= X(:)';
for j = 1:J
    len = ceil((length(A{j})+L-1)/2);
    for p = 1:(len-length(A{j+1}))
        [A,D] = post(j,A,D,H,G);
    end
end
A = A(2:end);
D = D(2:end);


function [S,A,D] = recursive(k,j,J,S,A,D,H,G)
if j > J
    return;
elseif mod(k,2) ~= 0
    [S,A,D] = update_dwt(j,S,A,D,H,G);
else
    [S,A,D] = recursive(k/2,j+1,J,S,A,D,H,G);
end

function [S,A,D] = update_dwt(j,S,A,D,H,G)
a = S(j,:)*H(end:-1:1);
d = S(j,:)*G(end:-1:1);
S(j+1,:) = push(S(j+1,:), a);
if A{j+1}(1) == 0
    A{j+1} = push(A{j+1}, a);
    D{j+1} = push(D{j+1}, d);
else
    A{j+1} = [A{j+1} a];
    D{j+1} = [D{j+1} d];
end

function Y = push(X,v)
Y = [X(length(v)+1:end) v];

function [A,D] = post(j,A,D,H,G)
k = length(A{j+1})+1;
i = 2*k-1;
d = max(0,i-length(A{j}));
m = 1:min(i,length(H)-d);
A{j+1}(k) = A{j}(i-d-m+1)*H(d+m);
D{j+1}(k) = A{j}(i-d-m+1)*G(d+m);