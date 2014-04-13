function X = rtidwt2(a,D,N,H,G)

count = 0;
J = size(D,1);
X = zeros(N,1);
L = length(H);
M = L/2;
C = 2*(M-1);
S = zeros(J,max(C,M+1));
D2 = zeros(J,max(C,M+1));
H = H(:);
G = G(:);

out = [];
S(end,end-C+1:end) = a(1:C);
D2(end,end-C+1:end) = D{end}(1:C);
for j = J:-1:1
    for p = C-M:-1:0
        if j == 1
            out = [out produce(p,j,S,D2,H,G)];
        else
            [S,D2,D] = update_idwt(p,j,S,D2,D,H,G);
        end
    end
end

for i = C+1:length(a)
    X = push(X, out);
    count = count + length(out);
    S(end,:) = push(S(end,:), a(i));
    D2(end,:) = push(D2(end,:), D{end}(i));
    [S,D2,D,out] = recursive2(0,J,S,D2,D,H,G);
end
X = push(X, out(1:N-count));


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