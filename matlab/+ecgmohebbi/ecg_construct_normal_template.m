function Result = ecg_construct_normal_template(Signal, Fs, R)
%   Controi um template de batida normal a partir da media das batidas
%   normais do ECG.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao dos picos de onda R
%
% Saída:
%   template de batida normal
%
N = length(Signal);
half = fix(0.4*Fs);
m = size(R,1);

Beats = zeros(m,half*2);
for i = 1:m
    left = max(1,R(i)-half);
    right = min(N,R(i)+half-1);
    Beat = Signal(left:right);
    samples = left+half-R(i)+(1:length(Beat));
    Beats(i,samples) = Beat';
end
Result = mean(Beats,1)';
