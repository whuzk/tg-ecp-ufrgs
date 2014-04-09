function [A,D] = rtdwt1(X,J,H,G)

L = length(H);
A = cell(J+1,1);
D = cell(J+1,1);
S = zeros(J+1,L);
H = H(:);
G = G(:);

len = 4+(L/2-1)*2^(J+1)-floor((L/2-1)/4)*2^J;
A{1} = [X(:)' zeros(1,len)];
for j = 1:J
    len = ceil(length(A{j})/2);
    A{j+1} = zeros(1,len);
    D{j+1} = zeros(1,len);
end

k = 1;
for i = 1:length(A{1})
    S(1,:) = push(S(1,:), A{1}(i));
    [S,A,D] = recursive(k,1,J,S,A,D,H,G);
    k = mod(k,2^J)+1;
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
A{j+1} = push(A{j+1}, a);
D{j+1} = push(D{j+1}, d);

function Y = push(X,v)
Y = [X(length(v)+1:end) v];