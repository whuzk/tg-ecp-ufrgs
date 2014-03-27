function Result = detect_qrs(Signal, Fs)
%   Detec�ao dos complexos QRS num sinal de ECG. Adaptado do metodo de
%   Tompkins.
%
% Entradas:
%   Signal - amplitudes normalizadas do sinal
%   Fs     - taxa de amostragem do sinal
%
% Sa�da:
%   localiza�ao dos picos de onda R
%

% filtering
[SignalF,SignalI,NewFs] = ecgfilter.tompkins_preprocess(Signal, Fs);

%{
[R,R2,THs,THn,THs2,THn2] = ...
    ecgfilter.tompkins_adapted(SignalF, SignalI, NewFs);
Result = round(R*Fs/NewFs);
%Result = ecgmath.neighbour_max(Signal,Result,10);
figure;
hold on, grid on;
plot([SignalF THs2 THn2]);
figure;
hold on, grid on;
plot([SignalI THs THn]);
plot(R, SignalI(R), 'kx');
plot(R2, SignalI(R2), 'ko');
ecgutilities.plot_signal_r(Signal, Result);
%}

% detection
R = ecgfilter.tompkins_production(SignalF, SignalI, NewFs);
Result = round(R*Fs/NewFs);