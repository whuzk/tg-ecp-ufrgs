function [A,D] = rtdwt2(X,J,H,G)

N = length(X);
L = length(H);
A = cell(J,1);
D = cell(J,1);
S = zeros(J+1,L);
len = N;
for j = 1:J
    len = floor(len/2);
    A{j} = zeros(1,len);
    D{j} = zeros(1,len);
end

k = 1;
for i = 1:length(X)
    S(1,:) = push(S(1,:), X(i));
    [S,A,D] = recursive(k,1,J,S,A,D,H,G);
    k = mod(k,2^J)+1;
end


function [S,A,D] = recursive(k,j,J,S,A,D,H,G)
if j > J
    return
elseif mod(k,2) == 0
    [S,A,D] = update_dwt(j,S,A,D,H,G);
    [S,A,D] = recursive(k/2,j+1,J,S,A,D,H,G);
end

function [S,A,D] = update_dwt(j,S,A,D,H,G)
a = S(j,:)*H(end:-1:1)';
d = S(j,:)*G(end:-1:1)';
S(j+1,:) = push(S(j+1,:), a);
A{j} = push(A{j}, a);
D{j} = push(D{j}, d);

function Y = push(X,v)
Y = [X(length(v)+1:end) v];