function [Filtered,Delay] = ecg_suppress_noise(Signal, Fs)
%   Remove o ruído de alta frequência de um sinal de ECG usando um filtro
%   de Butterworth.
%
% Entradas:
%   Signal - amplitudes do sinal original
%   Fs     - taxa de amostragem do sinal
%
% Saída:
%   Filtered - amplitudes do sinal filtrado
%   Delay - atraso médio da filtragem
%
N = 4;                          % ordem do filtro
Fc = 40;                        % frequencias de corte (em Hertz)
Wn = Fc * 2/Fs;                 % frequencias de corte normalizadas
[B,A] = butter(N, Wn);          % coeficientes do filtro
Filtered = filter(B, A, Signal);
Delay = fix(mean(grpdelay(B,A)));