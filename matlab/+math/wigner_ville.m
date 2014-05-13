function [Result,f] = wigner_ville(x, w)
% Obtém a transformada discreta de Wigner-Ville discreta

Lx = length(x);
N = length(w);
L1 = floor((N-1)/2);
L2 = ceil((N-1)/2);

x = math.hilbert_transform(x(:));
w = w(:);

Kernel = zeros(N,Lx);
for i = 1:Lx
    imin = min(i-1, Lx-i);
    tau = -min(imin,L1):min(imin,L2);
    index = rem(N+tau,N)+1;
    Kernel(index,i) = ...
        x(i+tau) .* conj(x(i-tau)) .* w(L1+1+tau) .* conj(w(L2+1-tau));
end;

Result = real(fft(Kernel));
f = 0.5*(0:N-1)'./N;