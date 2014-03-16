function [Filtered,Delay] = suppress_baseline(Signal, Fs)
%   Remove a linha de base flutuante de um sinal de ECG usando um filtro
%   wavelet adaptativo, de acordo com o metodo de Park.
%
% Entradas:
%   Signal - amplitudes do sinal original
%   Fs     - taxa de amostragem do sinal
%
% Sa�da:
%   Filtered - amplitudes do sinal filtrado�
%   Delay - atraso m�dio da filtragem
%
Filtered = Signal;
Delay = 0;