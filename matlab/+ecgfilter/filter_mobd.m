function [Result,delay] = filter_mobd(Signal,M)
N = length(Signal);
Result = zeros(N,1);
Signal = [ones(M-1,1); Signal(:)];
for i = M:N+M-1
    w = Signal(i-M+1:i);
    Result(i-M+1) = abs(prod(w));
end
delay = (M-1)/2;