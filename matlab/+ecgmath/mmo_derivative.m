function Result = mmo_derivative(Signal, B)
% Aplica a operaçao de derivaçao no sinal de acordo com um elemento
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
ke = M-1:-1:0;
kd = 0:M-1;
for i = M+delay:N
    e = min(Signal(i-ke)-B);
    d = max(Signal(i-kd)+B);
    f = Signal(i-delay);
    Result(i-delay) = (e+d-2*f)/delay;
end