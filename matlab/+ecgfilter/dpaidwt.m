function a = dpaidwt(a,D,Lx,H,G)

J = length(D);
L(J+1) = Lx;
for j = 1:J
    L(j) = length(D{j});
end

delay = floor(length(H)/2);
for j = 1:J
    a = upsample(a,2);
    d = upsample(D{j},2);
    a = wconv1(a(:),H) + wconv1(d(:),G);
    a = a(delay+(1:L(j+1)));
end