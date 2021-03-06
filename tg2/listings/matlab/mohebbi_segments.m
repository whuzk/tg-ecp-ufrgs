function Result = mohebbi_segments(Beats, J, Fs)

M = size(Beats,2);
L = round(0.08*Fs)*2;
Result = zeros(20,M);
for i = 1:M
    segment = Beats(J(i):J(i)+L-1,i);
    Result(:,i) = resample(segment,20,L);
end