function [F, D] = ecg_extract_gopalak_features(Signal, Fs, R, RR)
%   Extrai as caracteristicas do ECG pelo método de Mohebbi.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao das batidas
%   RR     - intervalos entre as batidas
%
% Saída:
%   F - matriz (M,l) com l = 0.8*Fs carateristicas para cada uma das M
%       batidas cardiacas extraidas do ECG
%   D - matriz (K,1) indicando quais M das K batidas foram detectadas
%
import ecggopalak.*;
global H Hs;

if isempty(R)
    error('no cardiac beats');
end

Hs = 2*fix(Fs*0.6)+1;
H = ecgmath.hermite(Hs, 50, 1.5);
D = true(size(R));
F = ecg_extract_hermite_coeff(Signal, Fs, R, min(RR,Hs));
