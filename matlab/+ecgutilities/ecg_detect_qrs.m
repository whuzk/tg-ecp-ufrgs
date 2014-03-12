function Result = ecg_detect_qrs(Signal, Fs)
%   Detecçao dos complexos QRS num sinal de ECG. Adaptado do metodo de
%   Tompkins.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%
% Saída:
%   localizaçao dos picos de onda R
%
N = length(Signal);
OriginalSignal = Signal;

%% filtering

% filtro passa-baixas (frequncia de corte aprox. 11Hz)
% low-pass (5T delay)
A1 = [1 -2 1];
B1([1 7 13]) = 1/36*[1 -2 1];
T1 = tf(B1,A1);

% high-pass (16T delay)
A2 = [1 -1];
B2([1 17 18 33]) = [-1/32 1 -1 1/32];
T2 = tf(B2,A2);

% derivative (2T delay)
A3 = 1;
B3 = 0.1*[2 1 0 -1 -2];
T3 = tf(B3,A3);

% cascade (23T delay)
T = T1*T2*T3;
Aa = fliplr(T.den{1});
Ba = T.num{1};

% integration (18T delay)
Ws = fix(0.15*Fs);
B4 = ones(1,Ws)/Ws;
T4 = tf(B4,1);

% smoothing (1T delay)
B5 = ones(1,3)/3;
T5 = tf(B5,1);

% cascade (19T delay)
T = T4*T5;
Ab = fliplr(T.den{1});
Bb = T.num{1};

% apply filters
Signal = filter(Ba, Aa, Signal);
Signal = filter(Bb, Ab, Signal.^2);
%ecgutilities.ecg_plot(Signal, 'integ');

%% detection
Stresh = 0;
Ntresh = 0;
Tresh1 = 0;
skip = false;
delay = 23+19;
M = fix(N/Fs*5);
R = zeros(M,1);
TTT = zeros(N,1);

k = 1;
i = 2;
while k <= M && i < N
    TTT(i) = Tresh1;
    if skip
        if Signal(i) < Tresh1
            skip = false;
        end
    else
        a = Signal(i-1);
        b = Signal(i);
        c = Signal(i+1);
        if (b-a) > 0 && (c-b) <= 0 && abs(c-2*b+a) > 1E-7
            PEAK = Signal(i);
            if PEAK > Tresh1
                Stresh = 0.25*PEAK + 0.75*Stresh;
                R(k) = i;
                k = k + 1;
                skip = true;
            else
                Ntresh = 0.25*PEAK + 0.75*Ntresh;
            end
            Tresh1 = Ntresh + 0.25*(Stresh-Ntresh);
        end
    end
    i = i + 1;
end
R = R(1:k-1);
Result = ecgmath.neighbour_max(OriginalSignal,R-delay,20);
%{
figure;
hold on, grid on;
plot([Signal TTT]);
plot(R, Signal(R), 'kx');
pause;
%}