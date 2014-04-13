function Result = rtclosing(Signal, B1, B2)

N = length(Signal);
M = length(B1);
if mod(M,2) == 0
    error('Structuring elements must have odd length.');
elseif length(B2) ~= M
    error('Structuring elements must have the same length.');
end

Result = zeros(N,1);
Signal = Signal(:);
Temp1 = zeros(M,1);
Temp2 = zeros(M,1);
B1 = wrev(B1(:));
B2 = B2(:);

delay = (M-1)/2;
k = 0;
for i = 1:N+2*delay
    k = mod(k,M)+1;
    if i <= N
        Temp1(k) = Signal(i);
    else
        Temp1(k) = 0;
    end
    if i > delay
        idx = [k+1:M 1:k];
        if i <= N+delay
            Temp2(k) = max(Temp1(idx)+B1);
        else
            Temp2(k) = 0;
        end
        if i > 2*delay
            Result(i-2*delay) = min(Temp2(idx)-B2);
        end
    end
end