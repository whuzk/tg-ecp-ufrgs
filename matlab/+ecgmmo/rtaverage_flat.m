function Result = rtaverage_flat(Signal, M1, M2)

N = length(Signal);
if mod(M1,2) == 0 || mod(M2,2) == 0
    error('Structuring elements must have odd length.');
end

Result = zeros(N,1);
Signal = Signal(:);
Temp = zeros(M1,1);
Temp12 = zeros(M1,1);
Temp13 = zeros(M2,1);
Temp14 = zeros(M2,1);
Temp22 = zeros(M1,1);
Temp23 = zeros(M2,1);
Temp24 = zeros(M2,1);

delay1 = (M1-1)/2;
delay2 = (M2-1)/2;
k1 = 0;
k2 = 0;
for i = 1:N+2*delay1+2*delay2
    k1 = mod(k1,M1)+1;
    k2 = mod(k2,M2)+1;
    if i <= N
        Temp(k1) = Signal(i);
    else
        Temp(k1) = 0;
    end
    if i > delay1
        idx1 = [k1+1:M1 1:k1];
        if i <= N+delay1
            Temp12(k1) = min(Temp(idx1));
            Temp22(k1) = max(Temp(idx1));
        else
            Temp12(k1) = 0;
            Temp22(k1) = 0;
        end
        if i > 2*delay1
            if i <= N+2*delay1
                Temp13(k2) = max(Temp12(idx1));
                Temp23(k2) = min(Temp22(idx1));
            else
                Temp13(k2) = 0;
                Temp23(k2) = 0;
            end
            if i > 2*delay1+delay2
                idx2 = [k2+1:M2 1:k2];
                if i <= N+2*delay1+delay2
                    Temp14(k2) = max(Temp13(idx2));
                    Temp24(k2) = min(Temp23(idx2));
                else
                    Temp14(k2) = 0;
                    Temp24(k2) = 0;
                end
                if i > 2*delay1+2*delay2
                    a = min(Temp14(idx2));
                    b = max(Temp24(idx2));
                    Result(i-2*delay1-2*delay2) = (a+b)/2;
                end
            end
        end
    end
end