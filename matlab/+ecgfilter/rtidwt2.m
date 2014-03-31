function X = rtidwt2(a,D,N,H,G)

J = size(D,1);
X = zeros(N,1);
L = length(H);
S = zeros(J,(L+2)/2);
D2 = zeros(J,(L+2)/2);

for i = 1:length(a)
    S(end,:) = push(S(end,:), a(i));
    D2(end,:) = push(D2(end,:), D{end}(i));
    [S,D2,D,out] = recursive(0,J,J,S,D2,D,H,G);
    X = push(X, out(:));
    for k = 1:2^(J-1)
        [S,D2,D,out] = recursive(k,1,J,S,D2,D,H,G);
        X = push(X, out(:));
    end
end


function [S,D2,D,out] = recursive(k,j,J,S,D2,D,H,G)
if k == 0
    if j == 1
        out = produce(0,j,S,D2,H,G);
    else
        [S,D2,D] = update_idwt((j<J),j,S,D2,D,H,G);
        [S,D2,D,out] = recursive(k,j-1,J,S,D2,D,H,G);
    end
elseif mod(k,2) ~= 0
    if j == 1
        out = produce(1-(k-1)/2,j,S,D2,H,G);
    else
        [S,D2,D] = update_idwt((k-1)/2,j,S,D2,D,H,G);
        out = [];
    end
else
    [S,D2,D,out] = recursive(k/2,mod(j,J-1)+1,J,S,D2,D,H,G);
end

function [S,D2,D] = update_idwt(p,j,S,D2,D,H,G)
out = produce(1-p,j,S,D2,H,G);
S(j-1,:) = push(S(j-1,:), out);
m = 1:length(H)/2;
D2(j-1,:) = push(D2(j-1,:), D{j-1}(m));
D{j-1}(m) = [];

function out = produce(p,j,S,D2,H,G)
m = 1:length(H)/2;
out(1) = S(j,end-m+p)*H(2*m-1)' + D2(j,end-m+p)*G(2*m-1)';
out(2) = S(j,end-m+p)*H(2*m)' + D2(j,end-m+p)*G(2*m)';

function X = push(X,v)
if isempty(v)
    return;
elseif isrow(X)
    X = [X(length(v)+1:end) v];
else
    X = [X(length(v)+1:end); v];
end
