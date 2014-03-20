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
[R,R1,Th] = ecgfilter.tompkins_decision(SignalB, SignalI, DelayI, NewFs);

% adjustment
R = round((R-DelayH)*Fs/NewFs);
Result = ecgmath.neighbour_max(Signal,R,10);

% plot
%
figure;
hold on, grid on;
plot([SignalI Th]);
plot(R1, SignalI(R1), 'kx');
%