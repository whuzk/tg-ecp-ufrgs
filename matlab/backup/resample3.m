function Result = resample3(Signal, NewFs, OldFs, h)

[P,Q] = rat(NewFs/OldFs, 1E-12);
Lx = length(Signal);
Ly = ceil(Lx*P/Q);
L = length(h);
delay = ceil(floor((L-1)/2)/Q);

Result = zeros(Ly,1);
temp = zeros(L,1);
q = 0;
l = 0;
n = -delay+1;
for i = 1:Lx+delay
    [out,temp,q,l,n] = res3(Signal,i,temp,h,q,l,n,P,Q,L,Lx,Ly);
    if ~isnan(out)
        Result(n) = out;
    end
end

function [out,temp,q,l,n] = res3(Signal,i,temp,h,q,l,n,P,Q,L,Lx,Ly)
out = NaN;
for p = 1:P
    if p == 1 && i < Lx
        x = Signal(i);
    else
        x = 0;
    end

    l = mod(l,L)+1;
    temp(l) = x;
    
    q = mod(q,Q)+1;
    if q == Q
        n = n + 1;
        if n > 0 && n < Ly
            out = h*[temp(l:-1:1); temp(end:-1:l+1)];
        end
    end
end