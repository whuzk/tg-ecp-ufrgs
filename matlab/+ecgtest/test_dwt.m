import ecgfilter.*

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
[A,D] = dpadwt(X,J,Lo_D,Hi_D);
[A1,D1] = rpadwt(X,J,Lo_D,Hi_D);
[A2,D2] = rpadwt2(X,J,Lo_D,Hi_D);

% Wavelet reconstruction
X1 = waverec(C,L,Lo_R,Hi_R);
X2 = dpaidwt(A{end},D,N,Lo_R,Hi_R);
X3 = rpaidwt(A1{end},D1,N,Lo_R,Hi_R);
X4 = dpaidwt(A2{end},D2,N,Lo_R,Hi_R);

%
% verification
%norm(X-X1)
%norm(X-X2)
%norm(X-X3)
%norm(X-X4)
%figure, plot([X X1]);
%figure, plot([X X2]);
%figure, plot([X X3]);
%figure, plot([X X4]);
%
figure, plot(X3);