% load mother wavelet
if ~exist('noisdopp','var')
    load noisdopp;
    %load leleccum;
    X = noisdopp(:);
    %X = leleccum(:);
    N = length(X);
    J = round(log2(N));
end

% load mother wavelet
wname = deblankl('db2');
[~,fname] = wavemngr('fields',wname,'type','file');
F = feval(fname,wname);

% normalize filter sum
W = F/sum(F);

% associated orthogonal filters
Lo_R = sqrt(2)*W;
Hi_R = qmf(Lo_R);
Hi_D = wrev(Hi_R);
Lo_D = wrev(Lo_R);

% wavelet decomposition
[C,L] = wavedec(X,J,Lo_D,Hi_D);
[A,D] = ecgmath.dpadwt(X,J,Lo_D,Hi_D);

% verification
norm(A{1}-C(1:L(1)));
end_i = L(1);
for i = 1:J
    norm(D{i}-C(end_i+(1:L(i+1))));
    end_i = end_i + L(i+1);
end

% Wavelet reconstruction
X1 = waverec(C,L,Lo_R,Hi_R);
X2 = ecgmath.dpaidwt(A{1},D,N,Lo_R,Hi_R);

% Verification
norm(X-X1)
norm(X-X2)
figure, plot([X X1]);
figure, plot([X X2]);
