function Result = mobd(Signal,M)
N = length(Signal);
Result = zeros(N,1);

a1 = floor((M-1)/2);
a2 = floor((M+1)/2);
Signal = [ones(a1,1); Signal; ones(a2,1);];

for i = M:N+M-1
    w = Signal(i-M+1:i);
    Result(i-M+1) = abs(prod(w));
end