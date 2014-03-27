function z = gpsmooth(y,wc,x)

% z=gpsmooth_3(y,wc,x) ... z is detrending (high pass) component
%
% INPUT
% y     - data
% wc    - normalised cutoff frequency ([0..1], where 1 corresponds to the ...
%         Nyquist frequency; 3dB attenuation, assuming the PSD is one at DC
% x     - time
%
% OUTPUT
% z     - filtered data
%
% AE 17/12/2008, 6/1/2009, 22/4/2009
%
% NOTE: This version of the code is different from the one reported
%         in the paper to run on MATLAB versions prior to 2006b.
%         Prior versions seem to have a bug in handling products of
%         scalars and sparse matrices.
%         Validated against 2006a ... 2008a
 
T=length(y);

% relation between smoothing parameter 3dB bandwidth
s2=(sqrt(2)-1)*(2*tan(wc*pi/2))^-4;
 
if nargin<4
    D2=spdiags(sqrt(s2)*ones(T-2,1)*[1 -2 1],[0:2],T-2,T);
else
    id=2:(T-1);
    idp1=id+1;
    idm1=id-1;
    V1=2./((x(idm1)-x(idp1)).*(x(id)-x(idp1)));
    V2=-2./((x(idm1)-x(id)).*(x(id)-x(idp1)));
    V3=2./((x(idm1)-x(id)).*(x(idm1)-x(idp1)));
    D2=spdiags(sqrt(s2)*[V1,V2,V3]/V1(1),[0:2],T-2,T);
end
 
% inverse kernel
invK=D2'*D2;
  
%evaluate the Cholesky factor R
R=chol(spdiags((1+sqrt(eps))*ones(T,1)*1,0,T,T)+invK);

% filter the data
z=R\(R'\y); 
