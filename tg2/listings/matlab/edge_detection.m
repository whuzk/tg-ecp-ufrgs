function pos = edge_detection(data, istart, iend, L, thr, default)

M = floor(L/2);
pos = default;
for i = istart+M:iend-M
    if all(data(i-M:i+M) < thr)
        pos = i;
        break;
    end
end