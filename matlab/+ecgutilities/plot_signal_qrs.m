function plot_signal_qrs(Filt,R1,R2,TH1,TH2,RR)

figure;
hold on, grid on;
plot([Filt TH1 TH2]);
plot(R1, Filt(R1), 'kx');
plot(R2, Filt(R2), 'ko');
figure, plot(RR);