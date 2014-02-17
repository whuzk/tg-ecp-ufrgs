function [A,I,J] = ecg_extract_st_deviation(Signal, Fs, R, RR)
%   Obtém os desvios de segmento ST de acordo com o método proposto por
%   Pang e com o método de analise tempo-frequencia proposto por Rocha.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localização das ondas R
%   RR     - intervalos entre as ondas R
%
% Saída:
%   A - desvios com base nos pontos de medida do método de Pang
%   B - desvios com base nos pontos de medida J e isoelétricos
%   J - localizaçao dos pontos J
%
import ecgmath.*;

N = length(Signal);
m = size(R,1);

half = fix(RR/2);
B1 = max(1,R-half);
B2 = min(N,R+half);

%% segmento ST baseado no pico da onda
HeartRate = 60*Fs ./ RR;
Limits = [-Inf 100 110 120 +Inf];
Offsets = fix([120 112 104 100]*1E-3 * Fs);

% obtem os pontos de medida
A = zeros(m,1);
for i = 1:length(Offsets)
    index = Limits(i) <= HeartRate & HeartRate < Limits(i+1);
    A(index) = min(R(index) + Offsets(i), N);
end

%% segmento ST baseado em analise tempo-frequencia
I = zeros(m,1);
J = zeros(m,1);
for i = 1:m
    Fres = 2^ceil(log2(0.1*Fs));
    Fmax = ceil(Fres*0.2);
    x = B1(i):B2(i);
    WVD = pseudo_wigner_ville(Signal(x), Fres);
    Sum = sum(WVD(1:Fmax,:),1) / Fmax;
    %Sum = sum(abs(WVD(1:Fmax,:)),1) / Fmax;
    %utilities.ecg_plot(Sum, 'Sum of Wigner-Ville transforms');
    
    % faixas de tempo
    t1 = (B1(i):R(i))-B1(i)+1;
    t2 = (R(i):B2(i))-B1(i)+1;
    
    % minimos locais
    [~,X1] = findpeaks(-Sum(t1));
    [~,X2] = findpeaks(-Sum(t2));
    X1 = X1 + B1(i) - 1;
    X2 = X2 + R(i) - 1;
    
    % ponto isoeletrico
    if ~isempty(X1)
        I(i) = X1(end);
    else
        I(i) = R(i) - fix(half(i)/6);
    end
    
    % ponto J
    if ~isempty(X2)
        J(i) = X2(1);
    else
        J(i) = R(i) + fix(half(i)/4);
    end
end
