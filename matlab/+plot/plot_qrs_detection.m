function plot_qrs_detection(sigI,qrs1,qrs2,th1,th2,rr1,rr2)
% Funçao para plotar a evoluçao do algoritmo de detecçao de batimentos

figure;
hold on, grid on;
plot([sigI th1 th2]);
plot(qrs1, sigI(qrs1), 'kx');
plot(qrs2, sigI(qrs2), 'ko');
title('Signal thresholds');

figure, plot([rr1 rr2]);
title('RR intervals');