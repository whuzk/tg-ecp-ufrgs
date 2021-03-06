function [S1,S2] = rocha_segments(Beats, B, J, E)

M = size(Beats,2);
S1 = zeros(64,M);
S2 = zeros(64,M);
for i = 1:M
    if (B(i) <= J(i))
        segment1 = Beats(B(i):J(i),i);
        S1(:,i) = resample(segment1,64,J(i)-B(i)+1);
    end
    if (J(i) <= E(i))
        segment2 = Beats(J(i):E(i),i);
        S2(:,i) = resample(segment2,64,E(i)-J(i)+1);
    end
end