function [F, D] = ecg_extract_rocha_features(Signal, Fs, R, RR)
%   Extrai as caracteristicas do ECG pelo método de Rocha.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao das batidas
%   RR     - intervalos entre as batidas
%
% Saída:
%   F - matriz (M,14) com 14 carateristicas para cada uma das M batidas
%       cardiacas extraidas do ECG
%   D - matriz (K,1) indicando quais M das K batidas foram detectadas
%
import ecgrocha.*;
import utilities.*;

if isempty(R)
    error('no cardiac beats');
end

D = true(size(R));
[A,I,J] = ecg_extract_st_deviation(Signal, Fs, R, RR);
[C1,C2] = ecg_hermite_coefficients(Signal, R, I, J, min(RR,Fs));

B1 = Signal(A) - Signal(I);
B2 = Signal(J) - Signal(I);
F = [B1 B2 C1 C2];

%ecg_plot_st_deviation(Signal, A, I, J);
%ecg_plot_hermite_expansion(Signal, R, I, J, min(RR,Fs), C1, C2);