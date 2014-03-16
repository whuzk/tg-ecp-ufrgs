function plot_signal(Signal, text)
% Funçao para plotar o grafico de um sinal de ECG

figure;
hold on;
grid on;
plot(Signal);
xlabel 'Sample';
ylabel 'Voltage (mV)';
title(text);
