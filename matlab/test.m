A = zeros(100,5000);
for i = 20:100
    pks = findpeaks(X,'minpeakdistance',i);
    M = 0;
    for j = 1:min(5000,length(pks))
        M = 0.9*M + 0.1*pks(j);
        A(i,j) = M;
    end
end
figure, surf(A);