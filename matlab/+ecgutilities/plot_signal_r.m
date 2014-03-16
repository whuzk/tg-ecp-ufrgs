function plot_signal_r(Signal, R)
% Funçao para plotar o grafico de um sinal de ECG com os picos de onda R

figure;
hold on;
grid on;
plot(Signal);
plot(R, Signal(R), 'kx');
xlabel 'Sample';
ylabel 'Voltage (mV)';
title('R peaks');
