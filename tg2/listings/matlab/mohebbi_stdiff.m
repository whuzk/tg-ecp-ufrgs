function C = mohebbi_stdiff(Sref, Stest)

[N,M] = size(Stest);
C = zeros(N,M);
for i = 1:M
    C(:,i) = Stest(:,i) - Sref;
end