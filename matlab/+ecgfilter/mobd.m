function Result = mobd(Signal,M)
N = length(Signal);
Result = zeros(N,1);

delay = floor((M-1)/2);
pad = ones(delay,1);
Signal = [pad; Signal; pad];

for i = M:N+M-1
    w = Signal(i-M+1:i);
    Result(i-M+1) = abs(prod(w));
end