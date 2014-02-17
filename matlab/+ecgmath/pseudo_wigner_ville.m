function Result = pseudo_wigner_ville(Signal, M)
%   Obtém a transformada de Wigner-Ville discreta com base numa janela
%   deslizante simétrica.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   M      - largura da janela deslizante
%
% Saída:
%   a transformada WVD do sinal
%
N = length(Signal);
L = round(M/2);
if (2^nextpow2(M) ~= M)
    warning('M should be a power of two.');
end;

Signal = hilbert(Signal);
Kernel = zeros(M,N);
for i = 1:N
    taumax = min([i-1, N-i, L-1]);
    tau = -taumax:taumax;
    index = rem(M+tau,M)+1;
    Kernel(index,i) = Signal(i+tau) .* conj(Signal(i-tau));
    
    if (L+1 <= i) && (i <= N-L)
        a = Signal(i+L) * conj(Signal(i-L));
        b = Signal(i-L) * conj(Signal(i+L));
        Kernel(L+1,i) = (a + b) / 2;
    end;
end;
Result = real(fft(Kernel));
