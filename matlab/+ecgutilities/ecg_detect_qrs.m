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
Stresh = 1E-4;
Ntresh = 1E-4;
achtung = false;
LASTPEAK = 0;
LASTDIFF = 0;
delay = 23+19;
M = fix(N/Fs*5);
R = zeros(M,1);
TTT = zeros(N,1);

k = 0;
i = 2;
while k <= M && i < N
    Tresh1 = Ntresh + 0.25*(Stresh-Ntresh);
    
    if Signal(i) > Tresh1
        if ~achtung
            achtung = true;
            if k == 0 || R(k) > 0
                k = k + 1;
            end
            %i, pause;
        end
    elseif achtung
        achtung = false;
        LASTPEAK = 0;
        %i, pause;
    end
    
    a = Signal(i-1);
    b = Signal(i);
    c = Signal(i+1);
    if (b-a) > 0 && (c-b) <= 0 && abs(c-2*b+a) > 1E-9
        PEAK = Signal(i);
        if achtung
            DIFF = PEAK - Tresh1;
            if DIFF > 0.2*LASTDIFF
                Stresh = 0.25*PEAK + 0.75*Stresh;
                if PEAK > LASTPEAK
                    LASTPEAK = PEAK;
                    LASTDIFF = DIFF;
                    R(k) = i;
                    %i, pause
                end
            end
        else
            Ntresh = 0.25*PEAK + 0.75*Ntresh;
        end
    end
    
    TTT(i) = Tresh1;
    i = i + 1;
end
R = R(1:k);
Result = ecgmath.neighbour_max(OriginalSignal,R-delay,10);
%{
figure;
hold on, grid on;
plot([Signal TTT]);
plot(R, Signal(R), 'kx');
%}