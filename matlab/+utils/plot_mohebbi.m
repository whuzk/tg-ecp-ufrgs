function plot_mohebbi(C, S)

M = size(C,2);
figure;
for i = 1:M
    plot(S(:,i)); hold on;
    plot(S(:,i) - C(:,i), '.k'); hold off;
    pause;
end