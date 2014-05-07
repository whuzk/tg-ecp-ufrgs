function plot_rpeaks(x,r)
% Funçao para plotar o grafico de um sinal de ECG com os picos de onda R

figure, plot(x);
hold on, grid on;
plot(r, x(r), 'xk');
xlabel 'Time (samples)';
ylabel 'Amplitude (adu/mV)';
title('R peak locations');
