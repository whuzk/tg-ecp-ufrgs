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
    wname = deblankl('db3');
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
%
% wavelet decomposition
[C,L] = wavedec(X,J,Lo_D,Hi_D);
[A,D] = dpadwt(X,J,Lo_D,Hi_D);
[A1,D1] = rpadwt(X,3,Lo_D,Hi_D);
[A2,D2] = rpadwt2(X,J,Lo_D,Hi_D);

% wavelet reconstruction
X1 = waverec(C,L,Lo_R,Hi_R);
X2 = dpaidwt(A{end},D,N,Lo_R,Hi_R);
X3 = rpaidwt(A1{end},D1,N,Lo_R,Hi_R);
X4 = dpaidwt(A2{end},D2,N,Lo_R,Hi_R);
%{
% de-noising
[XD,CXD,LXD] = wden(X,'sqtwolog','h','sln',J,wname);
XD2 = dpaidwt(A{end},denoise(D,J,N),N,Lo_R,Hi_R);
XD3 = rtdenoise(X,3,256,W);
figure, plot([XD2 XD3]);
%}
% verification
norm(X-X1)
norm(X-X2)
norm(X-X3)
norm(X-X4)
figure, plot([X X1]);
figure, plot([X X2]);
figure, plot([X X3]);
figure, plot([X X4]);
%
%% real-time
%{
[A3,D3] = rtdwt(X,3,Lo_D,Hi_D,Lo_R,Hi_R);
%
Signal1 = zeros(1,N);
Signal2 = zeros(1,N);
figure(1);
subplot(2,1,1); grid on;
title('Original signal');
lh1 = line((1:N)',Signal1);
xlim([1 N]);
subplot(2,1,2); grid on;
title('Reconstructed signal');
lh2 = line((1:N)',Signal2);
xlim([1 N]);
J = 3;
M = 256;
for i = 1:N
    Signal1 = [Signal1(2:end) X(i)];
    [A,D] = dpadwt(Signal1(end-M+1:end),J,Lo_D,Hi_D);
    XD = dpaidwt(A{end},denoise(D,J,M),M,Lo_R,Hi_R);
    %XD = dpaidwt(A{end},D,M,Lo_R,Hi_R);
    Signal2 = [Signal2(2:end) XD(end)];
    set(lh1, 'ydata', Signal1);
    set(lh2, 'ydata', Signal2);
    drawnow;
end
%norm(Signal2-Signal1)
%}