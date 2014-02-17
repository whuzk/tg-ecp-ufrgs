function Result = ecg_extract_artifact_beats(Signal, Fs, R)
%   Detecta batidas anormais presentes no sinal de ECG, com base numa
%   esimativa de presença dos pontos isoeletricos e pontos J.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao dos picos de onda R
%
% Saída:
%   indices das batidas classificadas como anormais
%
N = length(Signal);
m = size(R,1);
slope = 2.5/Fs;

Result = false(m,1);
for i = 1:m
    % ponto isoeletrico
    left = max(1,R(i)-fix(0.10*Fs));
    right = max(1,R(i)-fix(0.04*Fs));
    Window = Signal(left:right);
    ok1 = ~isempty(ecgmath.ecg_edge_detection(Window, slope));
    
    % ponto J
    left = min(N,R(i)+fix(0.02*Fs));
    right = min(N,R(i)+fix(0.12*Fs));
    Window = Signal(left:right);
    ok2 = ~isempty(ecgmath.ecg_edge_detection(Window, slope));
    
    % resultado
    Result(i) = ~ok1 && ~ok2;
end
