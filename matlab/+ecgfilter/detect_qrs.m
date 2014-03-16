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
[FiltSignal,Delay] = ecgfilter.tompkins_preprocess(Signal, Fs);
%ecgutilities.plot_signal(FiltSignal, 'integ');

% detection
%[~,R,~] = ecgfilter.tompkins_original(FiltSignal, Fs);
%ecgutilities.plot_signal_r(FiltSignal, R);
[R,Th] = ecgfilter.tompkins_adapted(FiltSignal, Fs);
%
figure;
hold on, grid on;
plot([FiltSignal Th]);
plot(R, FiltSignal(R), 'kx');
%

% adjustment
Result = ecgmath.neighbour_max(Signal,R-Delay,10);