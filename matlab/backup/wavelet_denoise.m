function D = wavelet_denoise(D,J,s)

%noise_level = median(abs(D{1}))/0.6745;
noise_level = std(D{1});
thr = s*noise_level;
for j = 1:J
    D{j} = D{j}.*(abs(D{j})>thr);
end