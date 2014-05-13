function plot_rocha(C, S1, S2)

M = size(C,2);
H1 = math.hermite_matrix(64,6,5);
H2 = math.hermite_matrix(64,6,8);
for i = 1:M
    S1hat = H1 * C(3:8,i);
    S2hat = H2 * C(9:14,i);
    figure(1);
    plot([S1(:,i); S2(:,i)]); hold on;
    plot([S1hat; S2hat], '.k'); hold off;
    stdev1 = [num2str(C(1,i)) ' (Pang)'];
    stdev2 = [num2str(C(2,i)) ' (Wigner-Ville)'];
    title(['Desvio ST: ' stdev1 ', ' stdev2]);
    pause;
end