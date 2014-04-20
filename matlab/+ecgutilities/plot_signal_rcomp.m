function plot_signal_rcomp(Signal,A,B,R)

scale2 = 0.1*max(Signal);
scale1 = 2*scale2;

figure;
hold on;
plot(Signal);
stem(R, A*scale1, 'g');
stem(R, B*scale2, 'r');