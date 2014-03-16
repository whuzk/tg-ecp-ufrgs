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

% filtering
[FiltSignal,Delay] = tompkins_filter(Signal, Fs);
%ecgutilities.ecg_plot(FiltSignal, 'integ');

% detection
R = sogari_qrs_detection(FiltSignal, Fs);
Result = ecgmath.neighbour_max(Signal,R-Delay,10);


function [Result,Delay] = tompkins_filter(Signal, Fs)
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

% smoothing (4T delay)
B5 = ones(1,9)/9;
T5 = tf(B5,1);

% smoothing (1T delay)
%B5 = ones(1,3)/3;
%T5 = tf(B5,1);

% cascade (22T delay)
T = T4*T5;
Ab = fliplr(T.den{1});
Bb = T.num{1};

% apply filters
Signal = filter(Ba, Aa, Signal);
Result = filter(Bb, Ab, Signal.^2);
Delay = 23 + 22;

function Result = sogari_qrs_detection(Signal, Fs)
%initializations
Sthresh = 0;
Nthresh = 0;
Dthresh = 0;
achtung = false;
LASTPEAK = 0;
N = length(Signal);
M = fix(N/Fs*5);
R = zeros(M,1);
TTT = zeros(N,1);
FOsize = 20;

% algorithm
k = 0;
i = 2;
while k <= M && i <= N-FOsize
    % get signal amplitude at current position
    P = Signal(i);
    
    % update global threshold according to the signal and noise thresholds
    Gthresh = Nthresh + 0.25*(Sthresh-Nthresh);
    
    % change state of the algorithm depending on whether the signal is
    % above the treshold or not
    if P > Gthresh
        if ~achtung
            % algorithm is entering the state of "active search" (when a
            % new QRS peak must be detected or refused)
            achtung = true;
            %i, pause;
        end
    elseif achtung
        % algorithm is leaving the state of active search
        achtung = false;
        % reset the last QRS peak (we know that the signal is in the range
        % of positive real values, so we can use zero as a special marker)
        LASTPEAK = 0;
        %i, pause;
    end
    
    % check if the current point is a positive peak or the verge of one
    a = Signal(i-1);
    c = Signal(i+1);
    if (P-a) > 0 && (c-P) <= 0 && abs(c-2*P+a) > 1E-9
        if ~achtung
            % peak has been considered as noise peak
            Nthresh = 0.125*P + 0.875*Nthresh;
            % give a decrease on the signal threshold, to quickly recover
            % from possible (and unusual) high peaks in the signal
            Sthresh = Sthresh - 0.025*(Sthresh-Nthresh);
            %i, pause
        elseif P-Gthresh < Dthresh
            % peak has been considered as artifact peak
            Dthresh = 0.5*Dthresh + 0.1*(P-Gthresh);
            %i, pause
        elseif P > Signal(i+FOsize)
            % peak has been considered as signal peak
            Sthresh = 0.125*P + 0.875*Sthresh;
            if P > LASTPEAK
                % peak is a new high for the current QRS
                Dthresh = 0.2*(P-Gthresh);
                if LASTPEAK == 0
                    % peak is the first for the current QRS
                    k = k + 1;
                end
                LASTPEAK = P;
                R(k) = i;
                %i, pause
            end
        end
    end
    
    TTT(i) = Gthresh;
    i = i + 1;
end
Result = R(1:k);
%
figure;
hold on, grid on;
plot([Signal TTT]);
plot(Result, Signal(Result), 'kx');
%