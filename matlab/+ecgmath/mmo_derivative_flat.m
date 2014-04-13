function Result = mmo_derivative_flat(Signal, s)
% Aplica a operaçao de derivaçao no sinal de acordo com um elemento
% estruturador nulo de tamanho 2*s+1

N = length(Signal);
Result = zeros(N,1);
Signal = Signal(:);

k = 0:2*s;
for i = 3*s+1:N
    w = Signal(i-k);
    f = Signal(i-s);
    Result(i-s) = (min(w)+max(w)-2*f)/s;
end