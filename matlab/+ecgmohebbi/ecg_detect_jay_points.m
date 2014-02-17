function Result = ecg_detect_jay_points(Signal, Fs, R)
%   Detecta os pontos J num sinal de ECG. Caso numa batida o ponto J nao
%   seja detectado, deixa o valor zero, indicando que nao esta presente.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao dos picos de onda R
%
% Saída:
%   lista com a localizaçao dos pontos J (ou valor zero para os que nao
%   foram detectados)
%
N = length(Signal);
m = size(R,1);
lf = fix(0.02*Fs);
l1 = fix(0.02*Fs);
l2 = fix(0.12*Fs);
B = ones(1,lf)/lf;
slope = 2.5/Fs;

Result = zeros(m,1);
for i = 1:m
    left = min(N,R(i)+l1);
    right = min(N,R(i)+l2);
    if left < right
        window = Signal(left:right);
        avg = filtfilt(B, 1, window);
        point = ecgmath.ecg_edge_detection(avg, slope);
        if ~isempty(point)
            Result(i) = point + left - 1;
        end
    end
end
