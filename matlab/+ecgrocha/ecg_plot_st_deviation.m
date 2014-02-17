function ecg_plot_st_deviation(Signal, A, I, J)

X = 1:length(Signal);

figure;
hold on, grid on;
plot(X, Signal, '-');
plot(A, Signal(A), 'ko');
plot(I, Signal(I), 'k^');
plot(J, Signal(J), 'ks');
title('ST-deviation points');
