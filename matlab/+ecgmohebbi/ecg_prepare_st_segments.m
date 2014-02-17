function Result = ecg_prepare_st_segments(Signal, J, TemplateST)
%   Obtem um conjunto de dados que caracteriza os segmentos ST de cada
%   batida e que possa ser usado como entrada de uma rede neural.
%
% Entradas:
%   Signal     - amplitudes normalizadas do sinal
%   J          - localizaçao dos pontos J
%   TemplateST - template do segmento ST
%
% Saída:
%   conjunto de dados dos segmentos ST
%
N = length(Signal);
m = size(J,1);
l1 = length(TemplateST);

% reduz o numero de amostras pela metade
model = (TemplateST(1:2:end) + TemplateST(2:2:end))/2;
l2 = length(model);

% faz a diferenca entre o segmento ST de cada batida e o modelo
Result = zeros(m,l2);
for i = 1:m
    ST = Signal(J(i):min(N,J(i)+l1-1));
    ST = (ST(1:2:end) + ST(2:2:end))/2;
    Result(i,:) = model - ST;
end
