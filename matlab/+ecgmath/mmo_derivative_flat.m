function Result = mmo_derivative_flat(Signal, s)
% Aplica a operaçao de derivaçao morfologica no sinal de acordo com um
% elemento estruturador constante igual a zero, de tamanho 2*s+1

Result = zeros(size(Signal));
Signal = [
    ones(s,1)*Signal(1)
    Signal(:)
    ones(s,1)*Signal(end)
];

for i = 1:length(Result)
    W = Signal(i:i+2*s);
    f = Signal(i+s);
    Result(i) = (min(W) + max(W) - 2*f) / s;
end
