function Result = mmo_derivative_flat(Signal, s)
% Aplica a operaçao de derivaçao no sinal de acordo com um elemento
% estruturador nulo de tamanho 2*s+1

B = ones(2*s+1,1);
Result = (imerode(Signal,B)+imdilate(Signal,B)-2*Signal)./s;
%{
N = length(Signal);
Result = zeros(N,1);
Signal = padarray(Signal(:),s,'symmetric','both');

k = 0:2*s;
for i = 2*s+1:N+2*s
    w = Signal(i-k);
    f = Signal(i-s);
    Result(i-2*s) = (min(w)+max(w)-2*f)/s;
end
%}