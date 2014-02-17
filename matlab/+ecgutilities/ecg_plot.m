function ecg_plot(Signal, text)
% Fun�ao para plotar o grafico de um ECG no plano Voltagem x Amostra

figure;
hold on;
grid on;
plot(Signal);
xlabel 'Sample';
ylabel 'Voltage (mV)';
title(text);
