function X = truertidwt1(a,D,N,H,G)

J = size(D,1);
X = zeros(N,1);
L = length(H);
S = zeros(J,L/2);
D2 = zeros(J,L/2);
H = H(:);
G = G(:);

count = 0;
Z = 2^(J-1);
i = 1;
while count < L-1
    if mod(i-Z,2^J) == 0
        [count,S,D2] = push_app(count,S,D2,a,D);
        if count > 1
            [S,D2,D] = update_idwt(mod(count,2),J-1,S,D2,D,H,G);
        end
    end
    i = i + 1;
end

k = 0;
for i = i:N+7*Z-1
    if mod(i-Z,2^J) == 0
        [count,S,D2] = push_app(count,S,D2,a,D);
    end
    if k < 2^J
        k = k + 1;
    else
        k = 1;
    end
    if i >= 6*Z
        [S,D2,D] = recursive(k,1,J,S,D2,D,H,G);
        out = produce(mod(k,2),1,S,D2,H,G);
        X = push(X, out);
    end
end


function [S,D2,D] = recursive(k,j,J,S,D2,D,H,G)
if j >= J
    return
elseif mod(k,2) ~= 0
    p = mod((k+1)/2+j,2);
    if mod(J,2) == 0 && j > 1
        p = ~p;
    end
    [S,D2,D] = update_idwt(p,j,S,D2,D,H,G);
else
    [S,D2,D] = recursive(k/2,j+1,J,S,D2,D,H,G);
end

function [S,D2,D] = update_idwt(p,j,S,D2,D,H,G)
out = produce(p,j+1,S,D2,H,G);
S(j,:) = push(S(j,:), out);
if ~isempty(D{j})
    D2(j,:) = push(D2(j,:), D{j}(1));
    D{j}(1) = [];
else
    D2(j,:) = push(D2(j,:), 0);
end

function out = produce(p,j,S,D2,H,G)
m = 1:length(H)/2;
out = S(j,end-m+1)*H(2*m-p) + D2(j,end-m+1)*G(2*m-p);

function [count,S,D2] = push_app(count,S,D2,a,D)
count = count + 1;
if count <= length(a)
    S(end,:) = push(S(end,:), a(count));
    D2(end,:) = push(D2(end,:), D{end}(count));
else
    S(end,:) = push(S(end,:), 0);
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