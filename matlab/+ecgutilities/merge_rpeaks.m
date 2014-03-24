function [A,B,R] = merge_rpeaks(R1, R2, Fs)

VDI1 = round(0.10*Fs);
VDI2 = round(0.20*Fs);

n = length(R1);
m = length(R2);
A = false(n+m,1);
B = false(n+m,1);
R = zeros(n+m,1);

i = 1;
j = 1;
k = 1;
while (i <= n) || (j <= m)
    while (i <= n) && (j <= m) && (R2(j) <= R1(i))
        if R1(i) - R2(j) <= VDI1
            A(k) = true;
            i = i + 1;
        end
        B(k) = true;
        R(k) = R2(j);
        k = k + 1;
        j = j + 1;
    end
    while (i <= n) && (j <= m) && (R1(i) <= R2(j))
        if R2(j) - R1(i) <= VDI2
            B(k) = true;
            j = j + 1;
        end
        A(k) = true;
        R(k) = R1(i);
        k = k + 1;
        i = i + 1;
    end
    while (i <= n) && (j > m)
        A(k) = true;
        R(k) = R1(i);
        k = k + 1;
        i = i + 1;
    end
    while (j <= m) && (i > n)
        B(k) = true;
        R(k) = R2(j);
        k = k + 1;
        j = j + 1;
    end
end

A(k:end) = [];
B(k:end) = [];
R(k:end) = [];