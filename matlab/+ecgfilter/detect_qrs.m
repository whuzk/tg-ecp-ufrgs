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
%ecgutilities.plot_signal(SignalI, 'integ');

% detection
%[R,Th] = ecgfilter.tompkins_adapted(SignalI, NewFs);
[R,R1,Th] = ecgfilter.tompkins_decision(SignalB, SignalI, DelayI, NewFs);
%
figure;
hold on, grid on;
plot([SignalI Th]);
plot(R1, SignalI(R1), 'kx');
%

% adjustment
%R = round((R-DelayH-DelayI)*Fs/NewFs);
R = round((R-DelayH)*Fs/NewFs);
Result = ecgmath.neighbour_max(Signal,R,10);