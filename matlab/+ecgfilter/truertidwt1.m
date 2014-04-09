function [X,delay] = truertidwt1(a,D,N,H,G)

J = size(D,1);
X = zeros(N,1);
L = length(H);
S2 = zeros(J+1,L/2);
D2 = zeros(J+1,L/2);
count = zeros(J+1,1);
H = H(:);
G = G(:);

i = 1;
while count(end) < L-1
    if mod(i,2^J) == 2^(J-1)
        [count(end),S2,D2] = push_app(count(end),S2,D2,a,D);
        save = count(1);
        [count,S2,D2,D] = recursive2(count,J,J,S2,D2,D,H,G);
        dif = count(1) - save;
        X = push(X, S2(1,end-dif+1:end));
    end
    i = i + 1;
end

k = 2^(J-1)+1;
delay = i+k-1+L+2;
for i = k:N-L+11
    if mod(i,2^J) == 0
        [count(end),S2,D2] = push_app(count(end),S2,D2,a,D);
    end
    [S2,D2,D] = recursive(k,1,J,S2,D2,D,H,G);
    [S2,D2,D] = update_idwt(mod(k,2),1,S2,D2,D,H,G);
    X = push(X, S2(1,end));
    k = mod(k,2^J)+1;
end


function [S2,D2,D] = recursive(k,j,J,S2,D2,D,H,G)
if j >= J
    return;
elseif mod(k,2) ~= 0
    p = mod((k+1)/2+j,2);
    [S2,D2,D] = update_idwt(p,j+1,S2,D2,D,H,G);
else
    [S2,D2,D] = recursive(k/2,j+1,J,S2,D2,D,H,G);
end

function [count,S2,D2,D] = recursive2(count,j,J,S2,D2,D,H,G)
if j == 0
    return;
elseif count(j+1) == length(H)/2 && count(j+1) < length(H)-1-(J-j)
    count(j) = count(j) + 1;
    [S2,D2,D] = update_idwt(0,j,S2,D2,D,H,G);
    [count,S2,D2,D] = recursive2(count,j-1,J,S2,D2,D,H,G);
elseif count(j+1) > length(H)/2
    count(j) = count(j) + 1;
    [S2,D2,D] = update_idwt(1,j,S2,D2,D,H,G);
    [count,S2,D2,D] = recursive2(count,j-1,J,S2,D2,D,H,G);
    if count(j+1) < length(H)-1-(J-j)
        count(j) = count(j) + 1;
        [S2,D2,D] = update_idwt(0,j,S2,D2,D,H,G);
        [count,S2,D2,D] = recursive2(count,j-1,J,S2,D2,D,H,G);
    end
end

function [S2,D2,D] = update_idwt(p,j,S2,D2,D,H,G)
out = produce(p,j+1,S2,D2,H,G);
S2(j,:) = push(S2(j,:), out);
if j > 1
    if ~isempty(D{j-1})
        D2(j,:) = push(D2(j,:), D{j-1}(1));
        D{j-1}(1) = [];
    else
        D2(j,:) = push(D2(j,:), 0);
    end
end

function out = produce(p,j,S2,D2,H,G)
m = 1:length(H)/2;
a = S2(j,end-m+1);
d = D2(j,end-m+1);
out = a*H(2*m-p) + d*G(2*m-p);

function [count,S2,D2] = push_app(count,S2,D2,a,D)
count = count + 1;
if count <= length(a)
    S2(end,:) = push(S2(end,:), a(count));
    D2(end,:) = push(D2(end,:), D{end}(count));
else
    S2(end,:) = push(S2(end,:), 0);
    D2(end,:) = push(D2(end,:), 0);
end

function X = push(X,v)
if isempty(v)
    return;
elseif isrow(X)
    X = [X(length(v)+1:end) v(:)'];
else
    X = [X(length(v)+1:end); v(:)];
end