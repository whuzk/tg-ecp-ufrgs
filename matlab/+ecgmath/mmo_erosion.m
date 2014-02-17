function Result = mmo_erosion(Signal, B)
% Aplica a opera�ao de erosao no sinal de acordo com um elemento
% estruturador B

N = length(Signal);
M = length(B);
l1 = floor((M-1)/2);
l2 = floor((M+1)/2);

Result = zeros(size(Signal));
Signal = [
    ones(l1,1)*Signal(1)
    Signal(:)
    ones(l2,1)*Signal(end)
];
B = B(:);

for i = 1:N
    W = Signal(i:i+M-1);
    Result(i) = min(W-B);
end
