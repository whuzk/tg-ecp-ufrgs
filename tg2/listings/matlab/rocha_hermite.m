function [C1,C2] = rocha_hermite(S1, S2)

[N,M] = size(S1);
C1 = zeros(6,M);
C2 = zeros(6,M);
H1 = math.hermite_matrix(N,6,5);
H2 = math.hermite_matrix(N,6,8);
H1PI = (H1'*H1)\H1';
H2PI = (H2'*H2)\H2';
for i = 1:M
    C1(:,i) = H1PI * S1(:,i);
    C2(:,i) = H2PI * S2(:,i);
end