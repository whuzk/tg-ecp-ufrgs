function Result = resample2(Signal, NewFs, OldFs, m, bta)

[p,q] = rat(NewFs/OldFs, 1E-12);
pqmax = max(p,q);
fc = 1/2/pqmax;
L = 2*m*pqmax+1;
f_r = [0 2*fc 2*fc 1];
a_r = [1 1 0 0];
h_r = firls(L-1, f_r, a_r);
h_r = p*h_r.*kaiser(L,bta)';
Lx = length(Signal);
delay = floor(L/2);

SignalUp(1:p:Lx*p) = Signal;
SignalFilt = conv2(SignalUp(:), h_r(:), 'full');
Result = SignalFilt(delay+1:q:end-delay);