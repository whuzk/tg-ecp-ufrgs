function C = rocha_stdev(Beats, I, J0, J1)

M = size(Beats,2);
C = zeros(2,M);
for i = 1:M
    C(1,i) = Beats(J0(i),i) - Beats(I(i),i);
    C(2,i) = Beats(J1(i),i) - Beats(I(i),i);
end