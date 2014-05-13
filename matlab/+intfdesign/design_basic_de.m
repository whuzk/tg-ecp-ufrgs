function [b,a,g,d] = design_basic_de(N,M)

b = nconv([1 -1],N);
c = nconv([1 1],M);
b = conv(b,c);
a = 1;
g = 2^max(N,M);
d = (N+M)/2;


function y = nconv(x,n)
y = 1;
for i = 1:n
    y = conv(y,x);
end