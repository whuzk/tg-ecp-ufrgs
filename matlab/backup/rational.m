function [p,q] = rational(x, tol)

nh = [1 0];    % [n(k) n(k-1)]
dh = [0 1];    % [d(k) d(k-1)]

k = 0;
save = x;
while (true)
    k = k + 1;
    d = round(x);
    x = x - d;
    
    save2 = nh(1);
    save3 = dh(1);
    nh(1) = nh(1) * d + nh(2);
    dh(1) = dh(1) * d + dh(2);
    nh(2) = save2;
    dh(2) = save3;
    
    dif = abs(nh(1) / dh(1) - save);
    if (x == 0 || dif <= tol)
        break;
    end
    x = 1 / x;
end
p = nh(1) / sign(dh(1));
q = abs(dh(1));