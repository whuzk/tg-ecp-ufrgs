function Result = rtdilation(Signal, B)

N = length(Signal);
M = length(B);
if mod(M,2) == 0
    error('Structuring element must have odd length.');
end

Result = zeros(N,1);
Signal = Signal(:);
Temp = zeros(M,1);
B = wrev(B(:));

delay = (M-1)/2;
k = 0;
for i = 1:N+delay
    k = mod(k,M)+1;
    if i <= N
        Temp(k) = Signal(i);
    else
        Temp(k) = 0;
    end
    if i > delay
        idx = [k+1:M 1:k];
        Result(i-delay) = max(Temp(idx)+B);
    end
end