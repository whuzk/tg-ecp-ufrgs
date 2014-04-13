function Result = rtopenclosing(Signal, B1, B2)

N = length(Signal);
M1 = length(B1);
M2 = length(B2);
if mod(M1,2) == 0 || mod(M2,2) == 0
    error('Structuring elements must have odd length.');
end

Result = zeros(N,1);
Signal = Signal(:);
Temp1 = zeros(M1,1);
Temp2 = zeros(M1,1);
Temp3 = zeros(M2,1);
Temp4 = zeros(M2,1);
B11 = B1(:);
B12 = wrev(B11);
B22 = B2(:);
B21 = wrev(B22);

delay1 = (M1-1)/2;
delay2 = (M2-1)/2;
k1 = 0;
k2 = 0;
for i = 1:N+2*delay1+2*delay2
    k1 = mod(k1,M1)+1;
    k2 = mod(k2,M2)+1;
    if i <= N
        Temp1(k1) = Signal(i);
    else
        Temp1(k1) = 0;
    end
    if i > delay1
        idx1 = [k1+1:M1 1:k1];
        if i <= N+delay1
            Temp2(k1) = min(Temp1(idx1)-B11);
        else
            Temp2(k1) = 0;
        end
        if i > 2*delay1
            if i <= N+2*delay1
                Temp3(k2) = max(Temp2(idx1)+B12);
            else
                Temp3(k2) = 0;
            end
            if i > 2*delay1+delay2
                idx2 = [k2+1:M2 1:k2];
                if i <= N+2*delay1+delay2
                    Temp4(k2) = max(Temp3(idx2)+B21);
                else
                    Temp4(k2) = 0;
                end
                if i > 2*delay1+2*delay2
                    Result(i-2*delay1-2*delay2) = min(Temp4(idx2)-B22);
                end
            end
        end
    end
end