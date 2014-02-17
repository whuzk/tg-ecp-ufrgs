function Result = ecg_filter(Signal, Fs)
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
N = 4;                                  % ordem do filtro
Fc = [0.5 40];                          % frequencias de corte (em Hertz)
Wn = Fc * 2 / Fs;                       % frequencias de corte normalizadas
[B,A] = butter(N, Wn);                  % coeficientes do filtro
Result = flipud(filter(B, A, Signal));  % aplica o filtro ao sinal
Result = flipud(filter(B, A, Result));  % (filtragem de fase zero)
[Peaks,~] = findpeaks(-Result);         % obtem os minimos da sinal
Result = Result - mean(-Peaks);         % subtrai pela media dos minimos
