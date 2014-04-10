function Y = rtdenoise3(X,J,M,wname)

%% intializations
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);

N = length(X);
L = length(Lo_D);
A = cell(J+1,1);
D = cell(J+1,1);
S = zeros(J+1,L);

Y = zeros(N,1);
S2 = zeros(J+1,L/2);
D2 = zeros(J+1,L/2);
count = zeros(J+1,1);
FL = zeros(M,1);

dwt_append_len = L*(2^J-1)-1;
idwt_prepend_len = L-min(7,floor(2^3));
out = zeros(1,idwt_prepend_len);
idwt_delay = (L-1+floor(J/5))*2^J;
total_len = N+dwt_append_len;
if J == 4
    idwt_delay = idwt_delay+4;
    total_len = total_len-4;
end

A{1} = [X(:)' zeros(1,dwt_append_len)];

siz = zeros(total_len,J+1);

%% DWT-IDWT
k1 = 1;
k2 = 1;
for i = 1:total_len
    % DWT
    S(1,:) = push(S(1,:), A{1}(i));
    [S,A,D] = recursive0(k1,1,J,S,A,D,Lo_D,Hi_D);
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
            [count,S2,D2,D,FL] = recursive2(count,J,J,S2,D2,D,Lo_R,Hi_R,FL);
            dif = count(1) - save;
            if dif > 0
                out(save+(1:dif)) = S2(1,end-dif+1:end);
            end
        end
        if i > idwt_delay-idwt_prepend_len
            m = i-idwt_delay+idwt_prepend_len;
            Y = push(Y, out(m));
        end
    else
        [S2,D2,D,FL] = recursive1(k2,1,J,S2,D2,D,Lo_R,Hi_R,FL);
        [S2,D2,D,FL] = update_idwt(mod(k2,2),1,S2,D2,D,Lo_R,Hi_R,FL);
        Y = push(Y, S2(1,end));
        k2 = mod(k2,2^J)+1;
    end
    
    for j = 1:J+1
        siz(i,j) = length(D{j});
    end
end

%figure, plot([siz(:,2) siz(:,3) siz(:,4)]);

function [S,A,D] = recursive0(k,j,J,S,A,D,H,G)
if j > J
    return;
elseif mod(k,2) ~= 0
    [S,A,D] = update_dwt(j,S,A,D,H,G);
else
    [S,A,D] = recursive0(k/2,j+1,J,S,A,D,H,G);
end

function [S,A,D] = update_dwt(j,S,A,D,H,G)
a = S(j,:)*H(end:-1:1)';
d = S(j,:)*G(end:-1:1)';
S(j+1,:) = push(S(j+1,:), a);
A{j+1} = [A{j+1} a];
D{j+1} = [D{j+1} d];

function [S2,D2,D,FL] = recursive1(k,j,J,S2,D2,D,H,G,FL)
if j >= J
    if j == 5 && k == 2
        [S2,D2,D,FL] = update_idwt(1,j,S2,D2,D,H,G,FL);
    end
    return;
elseif mod(k,2) ~= 0
    if j ~= 4 || k ~= 3
        p = mod((k+2*(j==3)-1)/2,2);
        [S2,D2,D,FL] = update_idwt(p,j+1,S2,D2,D,H,G,FL);
    end
else
    [S2,D2,D,FL] = recursive1(k/2,j+1,J,S2,D2,D,H,G,FL);
end

function [count,S2,D2,D,FL] = recursive2(count,j,J,S2,D2,D,H,G,FL)
if j == 0
    return;
elseif count(j+1) == length(H)/2
    if check(j-1,count(j),length(H))
        count(j) = count(j) + 1;
        [S2,D2,D,FL] = update_idwt(0,j,S2,D2,D,H,G,FL);
        [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
    end
elseif count(j+1) > length(H)/2
    if check(j-1,count(j),length(H))
        count(j) = count(j) + 1;
        [S2,D2,D,FL] = update_idwt(1,j,S2,D2,D,H,G,FL);
        [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
    end
    if check(j-1,count(j),length(H))
        count(j) = count(j) + 1;
        [S2,D2,D,FL] = update_idwt(0,j,S2,D2,D,H,G,FL);
        [count,S2,D2,D,FL] = recursive2(count,j-1,J,S2,D2,D,H,G,FL);
    end
end

function ok = check(j,count,L)
ok = (count < L-min(7,floor(2^(3-j))));

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
out = a*H(2*m-p)' + d.*(abs(d)>Thresh)*G(2*m-p)';

function [S2,D2,A,D] = push_app(S2,D2,A,D)
if ~isempty(A{end})
    S2(end,:) = push(S2(end,:), A{end}(1));
    D2(end,:) = push(D2(end,:), D{end}(1));
    A{end}(1) = [];
    D{end}(1) = [];
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