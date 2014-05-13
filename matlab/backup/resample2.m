function Result = resample2(Signal, NewFs, OldFs, h)

[p,q] = rat(NewFs/OldFs, 1E-12);
Lx = length(Signal);

SignalUp = zeros(Lx*p,1);
SignalUp(1:p:Lx*p) = Signal;
SignalFilt = conv2(SignalUp, h(:), 'same');
Result = SignalFilt(1:q:end);