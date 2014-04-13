function a = dpaidwt1(a,D,N,H,G)

J = length(D);
L = zeros(J,1);
L(1) = N;
for j = 2:J
    L(j) = length(D{j-1});
end

delay = length(H);
for j = J:-1:1
    a = dyadup(a(:),1,0);
    d = dyadup(D{j}(:),1,0);
    a = wconv1(a,H) + wconv1(d,G);
    a = a(delay+(1:L(j)));
end