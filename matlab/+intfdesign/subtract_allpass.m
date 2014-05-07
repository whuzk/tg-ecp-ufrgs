function [b,d] = subtract_allpass(b,a,g,d)

b1 = zeros(1,2*d+1);
d = ceil(d);
b1(d+1) = g;
b = conv(b1,a)-b;