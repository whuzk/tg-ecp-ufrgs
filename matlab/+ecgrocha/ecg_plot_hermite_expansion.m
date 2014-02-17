function ecg_plot_hermite_expansion(Signal, R, I, J, L, C1, C2)

import ecgmath.*;

N = length(Signal);
B = min(N,R+fix(L/2));

Hc = size(C1,2);
M1 = max(J-I+1);
M2 = max(B-J+1);
H1 = ecgmath.hermite_matrix(M1, Hc, 5);
H2 = ecgmath.hermite_matrix(M2, Hc, 8);

figure;
hold on, grid on;
plot(1:length(Signal), Signal);
for i = 1:size(R,1)
    X = I(i):J(i);
    len = length(X);
    h = (1:len)+fix((M1-len)/2);
    Y = H1(h,:)*C1(i,:)';
    plot(X, Y, 'k.');
    
    X = J(i):B(i);
    len = length(X);
    h = (1:len)+fix((M2-len)/2);
    Y = H2(h,:)*C2(i,:)';
    plot(X, Y, 'r.');
end
title('Hermitian expansion');
