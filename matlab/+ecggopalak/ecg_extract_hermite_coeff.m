function Result = ecg_extract_hermite_coeff(Signal, Fs, R, L)
%   Obtem os coeficientes das primeiras l funçoes discretas de Hermite na 
%   expansao de Hermite para as batidas de um sinal de ECG.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%   R      - localizaçao dos picos das ondas R
%   L      - largura das batidas cardiacas
%
% Saída:
%   matriz (m,l) com os coefficientes de Hermite
%
global H Hs;

N = length(Signal);
m = size(R,1);
half = fix(L/2);
Result = zeros(m,50);
for i = 1:m
    left = max(1,R(i)-half(i));
    right = min(N,R(i)+half(i));
    Beat = Signal(left:right);
    [x,y] = utilities.get_baseline_points(Beat, Fs, R(i)-left+1);
    Y = spline(x,y,(1:length(Beat))');
    X = ensure_length(Beat-Y, Hs);
    Result(i,:) = X'*H;
    %{
    figure(1), cla;
    hold on, grid on;
    plot(X, '-b');
    plot(H*Result(i,:)', '.k');
    pause;
    %}
end


function Result = ensure_length(Signal, l)

D = l-length(Signal);
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    Result = Signal(1+a:end-b);
else
    Result = [zeros(a,1); Signal; zeros(b,1)];
end
