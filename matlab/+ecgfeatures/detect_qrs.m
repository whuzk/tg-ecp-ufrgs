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
[SignalF,SignalI] = ecgfilter.tompkins_filter(Signal, Fs);
Result = ecgfeatures.tompkins_production(SignalF, SignalI, Fs);

%{
[R,R2,THs,THn,THs2,THn2] = ...
    ecgfeatures.tompkins_adapted(SignalF, SignalI, Fs);
Result = R;
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