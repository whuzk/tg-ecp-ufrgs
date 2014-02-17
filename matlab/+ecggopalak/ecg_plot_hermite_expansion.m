function ecg_plot_hermite_expansion(Signal, C, R, L)

N = length(Signal);
[m,l] = size(C);

half = fix(L/2);
Max = max(L);
h = fix(Max/2);
Max = 2*h+1;
j = -h:h;
u = ecgmath.hermite(Max, l, 1.5);

figure;
hold on, grid on;
plot(1:N, Signal, 'b-');
for i = 1:m
    k = R(i) + j;
    X = u*C(i,:)';
    a = (1 <= k & k <= N & -half(i) <= j & j <= half(i));
    plot(k(a), X(a), 'k.');
end
title('Hermite expansion');
