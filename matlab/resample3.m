function Result = resample3(Signal, NewFs, OldFs, m, bta)

[p,q] = rat(NewFs/OldFs, 1E-12);
pqmax = max(p,q);
fc = 1/2/pqmax;
L = 2*m*pqmax+1;
f_r = [0 2*fc 2*fc 1];
a_r = [1 1 0 0];
h_r = firls(L-1, f_r, a_r);
h_r = p*h_r.*kaiser(L,bta)';
Lx = length(Signal);
Ly = ceil(Lx*p/q);
delay = floor(L/2/q);

Result = zeros(Ly,1);
temp = zeros(1,L);
c = 0;
t = 0;
k = 0;
for i = -delay:Ly
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
        Result(i) = temp(t+1:end)*h_r(end:-1:t+1)' + temp(1:t)*h_r(t:-1:1)';
    end
end