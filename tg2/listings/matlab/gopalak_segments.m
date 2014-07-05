function Result = gopalak_segments(Beats, RR)

[N,M] = size(Beats);
RR = min(N,RR);
Result = zeros(250,M);
for i = 1:M
    off = floor((N-RR(i))/2);
    segment = Beats(off+1:off+RR(i),i);
    Result(:,i) = resample(segment,250,RR(i));
end