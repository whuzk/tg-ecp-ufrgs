function D = block_denoise(D,J,M)

for i = M:length(D{1})
    w = D{1}(max(1,i-M+1):i);
    thr = 4*std(w);
    k = i;
    j = 1;
    while mod(k,2) == 0 && j < J
        if abs(D{j}(k)) < thr
            D{j}(k) = 0;
        end
        k = k/2;
        j = j+1;
    end
    if abs(D{j}(k)) < thr
        D{j}(k) = 0;
    end
end