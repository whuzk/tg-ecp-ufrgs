function Result = resample3(Signal, NewFs, OldFs, h)

[p,q] = rat(NewFs/OldFs, 1E-12);
Lx = length(Signal);
Ly = ceil(Lx*p/q);
L = length(h);
delay = ceil(floor((L-1)/2)/q);

Result = zeros(Ly,1);
temp = zeros(L,1);
c = 0;
t = 0;
k = 0;
for i = -delay+2:Ly
    for j = 1:q
        t = mod(t,L)+1;
        c = mod(c,p)+1;
        if c == 1 && k < Lx
            k = k + 1;
            temp(t) = Signal(k);
        else
            temp(t) = 0;
        end
    end
    if i > 0
        Result(i) = h*[temp(t:-1:1); temp(end:-1:t+1)];
    end
end