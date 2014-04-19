function plot_signal_qrs(Signal,R1,R2,TH1,TH2,RR)

figure;
hold on, grid on;
plot([Signal TH1 TH2]);
plot(R1, Signal(R1), 'kx');
plot(R2, Signal(R2), 'ko');
figure, plot(RR);