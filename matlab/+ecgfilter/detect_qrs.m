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
[SignalB,SignalI,NewFs,delay] = ecgfilter.tompkins_preprocess(Signal, Fs);
%ecgutilities.plot_signal(FiltSignal, 'integ');

% detection
%[R,Th] = ecgfilter.tompkins_adapted(SignalI, NewFs);
[~,R] = ecgfilter.tompkins_decision(SignalB, SignalI, NewFs);
%{
figure;
hold on, grid on;
plot([SignalI Th]);
plot(R, SignalI(R), 'kx');
%}

% adjustment
%R = round((R-delay)*Fs/NewFs);
R = round((R-22)*Fs/NewFs);
Result = ecgmath.neighbour_max(Signal,R,10);