function Y = cyclicrtdenoise3(X,M,wname)

%% intializations
J = 3;

[~,fname] = wavemngr('fields',wname,'type','file');
F = feval(fname,wname);
F = F(:)/sum(F);

Lo_R = sqrt(2)*F;
Hi_R = qmf(Lo_R);
Hi_D = wrev(Hi_R);
Lo_D = wrev(Lo_R);

N = length(X);
L = length(F);
A = cell(J+1,1);
D = cell(J+1,1);
S = cell(J+1,1);
S2 = cell(J+1,1);
D2 = cell(J+1,1);
for j = 1:J+1
    S{j} = zeros(1,L);
    S2{j} = zeros(1,L/2);
    D2{j} = zeros(1,L/2);
end

Y = zeros(N,1);
FL = NaN(M,1);
count = zeros(J+1,1);
idx0 = min(1,2*(L/2-J));
idx1 = ones(J+1,1);
idx2 = ones(J+1,1);
idx3 = 1;

%dwt_append_len = 4+(L/2-1)*2^(J+1)-floor((L/2-1)/4)*2^J;
dwt_append_len = L*(2^J-1)-1;
total_len = N+dwt_append_len;
idwt_prepend_len = max(0,L-2^J+1);
out = zeros(1,idwt_prepend_len);
idwt_delay = (L-1)*2^J;

A{1} = [X(:)' zeros(1,dwt_append_len)];
A{2} = zeros(1,3*L+2);
A{3} = zeros(1,L+2);
A{4} = 0;
D{1} = [];
D{2} = zeros(1,3*L+2);
D{3} = zeros(1,L+2);
D{4} = 0;

%% DWT-IDWT
k1 = 1;
k2 = 2^(J-1)+1;

tic;
for i = 1:total_len
    % DWT
    S{1} = push(S{1}, A{1}(i));
    [S,A,D,idx1] = recursive0(k1,1,J,S,A,D,Lo_D,Hi_D,idx1);
    k1 = mod(k1,2^J)+1;
    
    % IDWT
    new_app = (mod(i,2^J) == 2^(J-1));
    if new_app
        [S2,D2,A,D] = push_app(S2,D2,A,D);
        count(end) = count(end) + 1;
    end
    if i <= idwt_delay
        if new_app
            save = count(1);
            [count,S2,D2,D,FL,idx2,idx3] = recursive2(count,J,J,S2,D2,D,Lo_R,Hi_R,FL,idx2,idx3);
            dif = count(1) - save;
            out(save+(1:dif)) = S2{1}(end-dif+1:end);
        end
        if i > idwt_delay-idwt_prepend_len
            m = i-idwt_delay+idwt_prepend_len;
            Y(idx0) = out(m);
            idx0 = idx0 + 1;
        end
    else
        [S2,D2,D,FL,idx2,idx3] = recursive1(k2,1,J,S2,D2,D,Lo_R,Hi_R,FL,idx2,idx3);
        [S2,D2,D,FL,idx2,idx3] = update_idwt(mod(k2,2),1,S2,D2,D,Lo_R,Hi_R,FL,idx2,idx3);
        if idx0 > 0
            Y(idx0) = S2{1}(end);
        end
        idx0 = idx0 + 1;
    end
    k2 = mod(k2,2^J)+1;
end
avg = toc/total_len;


function [S,A,D,idx] = recursive0(k,j,J,S,A,D,H,G,idx)
if j > J
    return;
elseif mod(k,2) ~= 0
    [S,A,D,idx] = update_dwt(j,S,A,D,H,G,idx);
else
    [S,A,D,idx] = recursive0(k/2,j+1,J,S,A,D,H,G,idx);
end

function [S,A,D,idx] = update_dwt(j,S,A,D,H,G,idx)
a = S{j}*H(end:-1:1);
d = S{j}*G(end:-1:1);
S{j+1} = push(S{j+1}, a);

k = idx(j+1);
A{j+1}(k) = a;
D{j+1}(k) = d;
idx(j+1) = mod(k,length(D{j+1}))+1;

function [S2,D2,D,FL,idx2,idx3] = recursive1(k,j,J,S2,D2,D,H,G,FL,idx2,idx3)
if j >= J
    return;
elseif mod(k,2) ~= 0
    p = mod((k+1)/2+j,2);
    [S2,D2,D,FL,idx2,idx3] = update_idwt(p,j+1,S2,D2,D,H,G,FL,idx2,idx3);
else
    [S2,D2,D,FL,idx2,idx3] = recursive1(k/2,j+1,J,S2,D2,D,H,G,FL,idx2,idx3);
end

function [count,S2,D2,D,FL,idx2,idx3] = recursive2(count,j,J,S2,D2,D,H,G,FL,idx2,idx3)
if j == 0
    return;
elseif count(j+1) == length(H)/2 && count(j+1) < length(H)-1-(J-j)
    count(j) = count(j) + 1;
    [S2,D2,D,FL,idx2,idx3] = update_idwt(0,j,S2,D2,D,H,G,FL,idx2,idx3);
    [count,S2,D2,D,FL,idx2,idx3] = recursive2(count,j-1,J,S2,D2,D,H,G,FL,idx2,idx3);
elseif count(j+1) > length(H)/2
    count(j) = count(j) + 1;
    [S2,D2,D,FL,idx2,idx3] = update_idwt(1,j,S2,D2,D,H,G,FL,idx2,idx3);
    [count,S2,D2,D,FL,idx2,idx3] = recursive2(count,j-1,J,S2,D2,D,H,G,FL,idx2,idx3);
    if count(j+1) < length(H)-1-(J-j)
        count(j) = count(j) + 1;
        [S2,D2,D,FL,idx2,idx3] = update_idwt(0,j,S2,D2,D,H,G,FL,idx2,idx3);
        [count,S2,D2,D,FL,idx2,idx3] = recursive2(count,j-1,J,S2,D2,D,H,G,FL,idx2,idx3);
    end
end

function [S2,D2,D,FL,idx2,idx3] = update_idwt(p,j,S2,D2,D,H,G,FL,idx2,idx3)
out = produce(p,j+1,S2,D2,H,G,FL);
S2{j} = push(S2{j}, out);
if j > 1
    k = idx2(j);
    D2{j} = push(D2{j}, D{j}(k));
    idx2(j) = mod(k,length(D{j}))+1;
    if j == 2
        FL(idx3) = D2{j}(end);
        idx3 = mod(idx3,length(FL))+1;
    end
end

function out = produce(p,j,S2,D2,H,G,FL)
if isnan(FL(end))
    ok = true(size(D2{j}));
else
    ok = abs(D2{j}) > 4*std(FL);
end
m = (length(H):-2:2)-p;
out = S2{j}*H(m) + D2{j}.*ok*G(m);

function [S2,D2,A,D] = push_app(S2,D2,A,D)
S2{end} = push(S2{end}, A{end});
D2{end} = push(D2{end}, D{end});

function X = push(X,v)
if isempty(v)
    return;
elseif isrow(X)
    X = [X(length(v)+1:end) v(:)'];
else
    X = [X(length(v)+1:end); v(:)];
end