function [A,B,R] = merge_rpeaks(R1, R2)

n = length(R1);
m = length(R2);
A = false(n+m,1);
B = false(n+m,1);
R = zeros(n+m,1);

i = 1;
j = 1;
k = 0;
while (i <= n) || (j <= m)
    while (i <= n) && (j <= m) && (R2(j) - R1(i) > 15)
        k = k + 1;
        A(k) = true;
        R(k) = R1(i);
        i = i + 1;
    end
    while (i <= n) && (j <= m) && (R1(i) - R2(j) > 15)
        k = k + 1;
        B(k) = true;
        R(k) = R2(j);
        j = j + 1;
    end
    if (i <= n) && (j <= m)
        if (abs(R2(j)-R1(i)) < 15)
            k = k + 1;
            A(k) = true;
            B(k) = true;
            R(k) = R2(j);
            i = i + 1;
            j = j + 1;
        elseif R1(i) < R2(j)
            i = i + 1;
        else
            j = j + 1;
        end
    elseif (i <= n)
        i = i + 1;
    elseif (j <= m)
        j = j + 1;
    end
end

A(k+1:end) = [];
B(k+1:end) = [];
R(k+1:end) = [];