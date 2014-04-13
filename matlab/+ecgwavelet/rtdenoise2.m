function Y = rtdenoise2(X,J,M,F)

%% intializations
Lo_R = sqrt(2)*F(:);
Hi_R = qmf(Lo_R);
Hi_D = wrev(Hi_R);
Lo_D = wrev(Lo_R);

N = length(X);
L = length(F);
A = cell(J+1,1);
D = cell(J+1,1);
S = zeros(J+1,L);

Y = zeros(N,1);
S2 = zeros(J+1,L/2);
D2 = zeros(J+1,L/2);
count = zeros(J+1,1);
FL = zeros(M,1);

append_len = 4+(L/2-1)*2^(J+1)-floor((L/2-1)/4)*2^J;
dwt_len = N+append_len;
idwt_len = N+L*(2^J-1)-1;

%% DWT
A{1} = [X(:)' zeros(1,append_len)];
for j = 1:J
    len = ceil(length(A{j})/2);
    A{j+1} = zeros(1,len);
    D{j+1} = zeros(1,len);
end

k = 1;
for i = 1:dwt_len
    S(1,:) = push(S(1,:), A{1}(i));
    [S,A,D] = recursive0(k,1,J,S,A,D,Lo_D,Hi_D);
    k = mod(k,2^J)+1;
end

%% IDWT + denoise
prepend_len = max(0,L-2^J+1);
out = zeros(1,prepend_len);
delay = (L-1)*2^J;
a = A{end};
for i = 1:delay
    if mod(i,2^J) == 2^(J-1)
        [count(end),S2,D2] = push_app(count(end),S2,D2,a,D);
        save = count(1);
        [count,S2,D2,D,FL] = recursive2(count,J,J,S2,D2,D,Lo_R,Hi_R,FL);
        dif = count(1) - save;
        out(save+(1:dif)) = S2(1,end-dif+1:end);
    end
    if i > delay-prepend_len
        m = i-delay+prepend_len;
        Y = push(Y, out(m));
    end
end

k = 2^(J-1)+1;
for i = 1:idwt_len-delay
    if mod(i,2^J) == 2^(J-1)
        [count(end),S2,D2] = push_app(count(end),S2,D2,a,D);
    end
    [S2,D2,D,FL] = recursive1(k,1,J,S2,D2,D,Lo_R,Hi_R,FL);
    [S2,D2,D,FL] = update_idwt(mod(k,2),1,S2,D2,D,Lo_R,Hi_R,FL);
    Y = push(Y, S2(1,end));
    k = mod(k,2^J)+1;
end


function [S,A,D] = recursive0(k,j,J,S,A,D,H,G)
if j > J
    return;
elseif mod(k,2) ~= 0
    [S,A,D] = update_dwt(j,S,A,D,H,G);
else
    [S,A,D] = recursive0(k/2,j+1,J,S,A,D,H,G);
end

function [S,A,D] = update_dwt(j,S,A,D,H,G)
a = S(j,:)*H(end:-1:1);
d = S(j,:)*G(end:-1:1);
S(j+1,:) = push(S(j+1,:), a);
A{j+1} = push(A{j+1}, a);
D{j+1} = push(D{j+1}, d);

function [S2,D2,D,FL] = recursive1(k,j,J,S2,D2,D,H,G,FL)
if j >= J
    return;
elseif mod(k,2) ~= 0
    p = mod((k+1)/2+j,2);
    [S2,D2,D,FL] = update_idwt(p,j+1,S2,D2,D,H,G,FL);
else
    [S2,D2,D,FL] = recursive1(k/2,j+1,J,S2,D2,D,H,G,FL);
end

function [count,S2,D2,D,FL] = recursive2(count,j,J,S2,D2,D,H,G,FL)
if j == 0
    return;
elseif count(j+1) == length(H)/2 && count(j+1) < length(H)-1-(J-j)
    count(j) = count(j) + 1;
    [S2,D2,D,FL] = update_idwt(0,j,S2,D2,D,H,G,FL);
    [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
elseif count(j+1) > length(H)/2
    count(j) = count(j) + 1;
    [S2,D2,D,FL] = update_idwt(1,j,S2,D2,D,H,G,FL);
    [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
    if count(j+1) < length(H)-1-(J-j)
        count(j) = count(j) + 1;
        [S2,D2,D,FL] = update_idwt(0,j,S2,D2,D,H,G,FL);
        [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
    end
end

function [S2,D2,D,FL] = update_idwt(p,j,S2,D2,D,H,G,FL)
out = produce(p,j+1,S2,D2,H,G,FL);
S2(j,:) = push(S2(j,:), out);
if j > 1
    if ~isempty(D{j})
        D2(j,:) = push(D2(j,:), D{j}(1));
        D{j}(1) = [];
    else
        D2(j,:) = push(D2(j,:), 0);
    end
    if j == 2
        FL = push(FL, D2(j,end));
    end
end

function out = produce(p,j,S2,D2,H,G,FL)
Thresh = 0;%4*std(FL);
m = 1:length(H)/2;
a = S2(j,end-m+1);
d = D2(j,end-m+1);
out = a*H(2*m-p) + d.*(abs(d)>Thresh)*G(2*m-p);

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