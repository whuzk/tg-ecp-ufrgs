function [S1,S2] = rocha_segments(Beats, B, J, E, Fs)

M = size(Beats,2);
S1 = zeros(64,M);
S2 = zeros(64,M);
for i = 1:M
    if (B(i) <= J(i))
        segment1 = Beats(B(i):J(i),i);
        sresamp1 = resample(segment1,64,J(i)-B(i)+1);
        S1(:,i) = sresamp1;
    end
    if (J(i) <= E(i))
        segment2 = Beats(J(i):E(i),i);
        sresamp2 = resample(segment2,64,E(i)-J(i)+1);
        S2(:,i) = sresamp2;
    end
    %figure(1), plot([segment1; segment2]);
    %figure(2), plot([sresamp1; sresamp2]);
    %pause;
end