function C = gopalak_hermite(Segments)

[N,M] = size(Segments);
C = zeros(50,M);
H = math.hermite(N,50,1.5)';
for i = 1:M
    C(:,i) = H * Segments(:,i);
end