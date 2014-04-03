function X = rtidwt2(a,D,N,H,G)

J = size(D,1);
X = zeros(N,1);
L = length(H);
S = zeros(J,(L+2)/2);
D2 = zeros(J,(L+2)/2);
H = H(:);
G = G(:);

S(end,end-1:end) = a(1:2);
D2(end,end-1:end) = D{end}(1:2);
[S,D2,D,out] = recursive1(0,J,S,D2,D,H,G);
for i = 3:length(a)
    X = push(X, out);
    S(end,:) = push(S(end,:), a(i));
    D2(end,:) = push(D2(end,:), D{end}(i));
    [S,D2,D,out] = recursive2(0,J,S,D2,D,H,G);
end
X = push(X, out(1:end-2));


function [S,D2,D,out] = recursive1(p,j,S,D2,D,H,G)
if j == 1
    out = produce(p,j,S,D2,H,G);
else
    [S,D2,D] = update_idwt(p,j,S,D2,D,H,G);
    [S,D2,D,out] = recursive1(p,j-1,S,D2,D,H,G);
end

function [S,D2,D,out] = recursive2(p,j,S,D2,D,H,G)
if j == 1
    out = produce(p,j,S,D2,H,G);
else
    [S,D2,D] = update_idwt(p,j,S,D2,D,H,G);
    [S,D2,D,out1] = recursive2(1,j-1,S,D2,D,H,G);
    [S,D2,D,out2] = recursive2(0,j-1,S,D2,D,H,G);
    out = [out1 out2];
end

function [S,D2,D] = update_idwt(p,j,S,D2,D,H,G)
out = produce(p,j,S,D2,H,G);
S(j-1,:) = push(S(j-1,:), out);
l = min(length(D{j-1}),2);
newD = [D{j-1}(1:l) zeros(1,2-l)];
D2(j-1,:) = push(D2(j-1,:), newD);
D{j-1}(1:l) = [];

function out = produce(p,j,S,D2,H,G)
m = 1:length(H)/2;
out(1) = S(j,end-m+1-p)*H(2*m-1) + D2(j,end-m+1-p)*G(2*m-1);
out(2) = S(j,end-m+1-p)*H(2*m) + D2(j,end-m+1-p)*G(2*m);

function X = push(X,v)
if isempty(v)
    return;
elseif isrow(X)
    X = [X(length(v)+1:end) v(:)'];
else
    X = [X(length(v)+1:end); v(:)];
end
