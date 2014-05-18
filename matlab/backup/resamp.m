function [y,h] = resamp(x,p,q)

Lx = length(x)
[p,q] = rat(p/q, 1e-12)

pqmax = max(p,q)
wc = pi/pqmax
M = 500
M2 = 250
%M = 2*10*pqmax;
%M2 = 10*pqmax;

n = (0:M)';
w = 0.54-0.46*cos(2*pi*n/M);
h = sin(wc*(n-M2))./(pi*(n-M2));
h(M2+1) = wc/pi;
h = p.*h.*w;
y = upfirdn(x,h,p,q);

delay = floor(M2/q)
Ly = ceil(Lx*p/q)
y = y(delay+(1:Ly));