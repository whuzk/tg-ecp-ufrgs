function [b,a,g,d] = design_basic_lp(N,m)

b = [1 zeros(1,m-1) -1];
a = [1 -1];
b = nconv(b,N);
a = nconv(a,N);
g = m^N;
d = N*(m-1)/2;


function y = nconv(x,n)
y = 1;
for i = 1:n
    y = conv(y,x);
end