import ecgfilter.*

% load mother wavelet
if ~exist('noisdopp','var')
    load noisdopp;
    X = noisdopp(:);
    %load leleccum;
    %X = leleccum(:);
    N = length(X);
    J = round(log2(N));
end

if ~exist('F','var')
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
end

%% analytic
% wavelet decomposition
[C,L] = wavedec(X,J,Lo_D,Hi_D);
[A1,D1] = dpadwt1(X,J,Lo_D,Hi_D);
[A2,D2] = dpadwt2(X,J,Lo_D,Hi_D);
[A3,D3] = rpadwt1(X,J,Lo_D,Hi_D);
[A4,D4] = rpadwt2(X,J,Lo_D,Hi_D);
[A5,D5] = rtdwt2(X,J,Lo_D,Hi_D);

% wavelet reconstruction
XR = waverec(C,L,Lo_R,Hi_R);
X1 = dpaidwt1(A1{end},D1,N,Lo_R,Hi_R);
X2 = dpaidwt2(A2{end},D2,N,Lo_R,Hi_R);
X3 = dpaidwt1(A3{end},D3,N,Lo_R,Hi_R);
X4 = dpaidwt2(A4{end},D4,N,Lo_R,Hi_R);
X5 = rtidwt2(A5{end},D5,N,Lo_R,Hi_R);

% verification
norm(X-XR)
norm(X-X1)
norm(X-X2)
norm(X-X3)
norm(X-X4)
norm(X-X5)
figure, plot([X XR]);
figure, plot([X X1]);
figure, plot([X X2]);
figure, plot([X X3]);
figure, plot([X X4]);
figure, plot([X X5]);
%{
% de-noising
XD = wden(X,'sqtwolog','h','sln',J,wname);
XD1 = dpaidwt1(A1{end},denoise(D1,J,N),N,Lo_R,Hi_R);
XD2 = dpaidwt2(A2{end},denoise(D2,J,N),N,Lo_R,Hi_R);
XD3 = dpaidwt1(A3{end},denoise(D3,J,N),N,Lo_R,Hi_R);
XD4 = dpaidwt2(A4{end},denoise(D4,J,N),N,Lo_R,Hi_R);
XD5 = rtidwt2(A5{end},denoise(D5,J,N),N,Lo_R,Hi_R);

norm(XD-XD1)
norm(XD-XD2)
norm(XD-XD3)
norm(XD-XD4)
norm(XD-XD5)
figure, plot([XD XD1]);
figure, plot([XD XD2]);
figure, plot([XD XD3]);
figure, plot([XD XD4]);
figure, plot([XD XD5]);
%}