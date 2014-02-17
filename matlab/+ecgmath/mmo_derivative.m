function Result = mmo_derivative(Signal, Gs)
% Aplica a operaçao de derivaçao morfologica no sinal de acordo com um
% elemento estruturador Gs (Gs deve ter largura 2*s+1)

N = length(Signal);
M = length(Gs);
l1 = floor((M-1)/2);
l2 = floor((M+1)/2);

Result = zeros(size(Signal));
Signal = [
    ones(l1,1)*Signal(1)
    Signal(:)
    ones(l2,1)*Signal(end)
];
Gs = Gs(:);
s = l1;

for i = 1:N
    W = Signal(i:i+M-1);
    a = min(W-Gs);
    b = max(W+Gs);
    f = Signal(i+s);
    Result(i) = (a + b - 2*f) / s;
end
