function Result = ecg_suppress_noise(Signal, Fs)
%   Faz a filtragem de um sinal de ECG usando um filtro de Butterworth e
%   elimina a linha de base.
%
% Entradas:
%   Signal - amplitudes do sinal original
%   Fs     - taxa de amostragem do sinal
%
% Saída:
%   amplitudes do sinal filtrado
%
N = 4;                          % ordem do filtro
Fc = 40;                        % frequencias de corte (em Hertz)
Wn = Fc * 2/Fs;                 % frequencias de corte normalizadas
[B,A] = butter(N, Wn);          % coeficientes do filtro
Result = filter(B, A, Signal);
