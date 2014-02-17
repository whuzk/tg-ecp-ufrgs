function [C1,C2] = ecg_hermite_coefficients(Signal, R, I, J, L)
%   Obtem os coeficientes das funçoes de hermite que melhor caracterizam os
%   complexos QRS e as ondas T de um ECG, de acordo com o metodo de Rocha.
%   Cada batida é dividida em dois segmentos, um contendo o complexo QRS e
%   outro contendo a onda T. A regressao é feita separadamente entre os
%   segmentos. Apenas os coeficientes das Hc primeiras funçoes de hermite
%   sao considerados.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localização dos picos de onda R
%   I      - localização dos pontos isoeletricos
%   J      - localizaçao dos pontos J
%   L      - larguras das batidas cardiacas
%
% Saída:
%   C1 - matriz (M,6) de coeficientes das seis primeiras funçoes de hermite
%        dos segmentos contendo o complexo QRS, para cada uma das M batidas
%   C2 - matriz (M,6) de coeficientes das seis primeiras funçoes de hermite
%        dos segmentos contendo a onda T, para cada uma das M batidas
%
Hc = 6;
A1 = 5;
A2 = 8;

N = length(Signal);
m = size(R,1);
C1 = zeros(m,Hc);
C2 = zeros(m,Hc);
B = min(N,R+fix(L/2));

M1 = max(J-I+1);
M2 = max(B-J+1);
H1 = ecgmath.hermite_matrix(M1, Hc, A1);
H2 = ecgmath.hermite_matrix(M2, Hc, A2);
for i = 1:m
    % segmento contendo o complexo QRS
    if J(i) - I(i) >= Hc
        X = I(i):J(i);
        Y = Signal(X);
        len = length(Y);
        h = (1:len) + fix((M1-len)/2);
        H = H1(h,:);
        C1(i,:) = (H'*H)\H' * Y;
    end
    
    % segmento contendo a onda T
    if B(i) - J(i) >= Hc
        X = J(i):B(i);
        Y = Signal(X);
        len = length(Y);
        h = (1:len) + fix((M2-len)/2);
        H = H2(h,:);
        C2(i,:) = (H'*H)\H' * Y;
    end
end
