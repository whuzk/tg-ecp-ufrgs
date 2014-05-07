function [b,a,g,d] = design_basic_bp(N,m,theta)

b = [1 zeros(1,m-1) -1];
a = [1 round(-2*cos(theta)) 1];
b = nconv(b,N);
a = nconv(a,N);
g = (m/2*abs(cos(m/2*theta))/sin(theta))^N;
d = N*(m/2-1);


function y = nconv(x,n)
y = 1;
for i = 1:n
    y = conv(y,x);
end