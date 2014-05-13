function plot_rpeak_comparison(data,a,b,r)
% Funçao para plotar o grafico de um sinal de ECG com os picos de onda R
% conhecidos e os picos detectados pelo algoritmo de detecçao

scale2 = 0.1*max(data);
scale1 = 2*scale2;

figure, plot(data);
hold on, grid on;
stem(r, a*scale1, 'g');
stem(r, b*scale2, 'r');
title('Annotated QRS (green) vs Detected QRS (red)');