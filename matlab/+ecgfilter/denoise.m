function D = denoise(D,J,N)

noise_level = std(D{1});
threshold_s = 4;

%noise_level = median(abs(D{1}))/0.6745;
%threshold_s = sqrt(2*log(N));

thr = threshold_s*noise_level;
for j = 1:J
    D{j} = D{j}.*(abs(D{j})>thr);
end