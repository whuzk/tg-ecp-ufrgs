function a = dpaidwt(a,D,Lx,H,G)

J = length(D);
for j = 1:J
    a = dyadup(a,0);
    d = dyadup(D{j},0);
    a = wconv1(a,H) + wconv1(d,G);
    if j < J
        a = wkeep1(a,length(D{j+1}));
    else
        a = wkeep1(a,Lx);
    end
end