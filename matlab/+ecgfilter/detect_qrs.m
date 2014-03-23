function Result = detect_qrs(Signal, Fs)
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
[SignalB,SignalI,DelayH,DelayI,NewFs] = ecgfilter.tompkins_preprocess(Signal, Fs);

% detection
%[R,R2,THs,THn,THs2,THn2] ...
% = ecgfilter.tompkins_adapted(SignalB, SignalI, DelayI, NewFs);
R = ecgfilter.tompkins_production(SignalB, SignalI, DelayI, NewFs);

% adjustment
Radj = round((R-DelayH-DelayI)*Fs/NewFs);
Result = ecgmath.neighbour_max(Signal,Radj,10);

%{
% plot
figure;
hold on, grid on;
plot([SignalB THs2 THn2]);
figure;
hold on, grid on;
plot([SignalI THs THn]);
plot(R, SignalI(R), 'kx');
plot(R2, SignalI(R2), 'ko');
ecgutilities.plot_signal_r(Signal, Result);
%}