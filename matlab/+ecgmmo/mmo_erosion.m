function Result = mmo_erosion(Signal, B)
% Aplica a operaçao de erosao no sinal de acordo com um elemento
% estruturador B

N = length(Signal);
M = length(B);
if mod(M,2) == 0
    error('Structuring element B must have odd length.');
end

Result = zeros(N,1);
Signal = Signal(:);
B = B(:);

delay = (M-1)/2;
pad = zeros(delay,1);
Signal = [pad; Signal; pad];

k = M-1:-1:0;
for i = M:N+M-1
    Result(i-M+1) = min(Signal(i-k)-B);
end