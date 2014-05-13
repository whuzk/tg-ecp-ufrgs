function plot_gopalak_features(C, S)

M = size(C,2);
H = math.hermite(250,50,1.5);
for i = 1:M
    Shat = H * C(:,i);
    figure(1);
    plot(S(:,i)); hold on;
    plot(Shat, '.k'); hold off;
    pause;
end