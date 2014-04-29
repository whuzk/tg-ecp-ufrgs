function Result = frame_beat(Beat, L)

D = L-length(Beat);
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    Result = Beat(1+a:end-b);
else
    Result = [zeros(a,1); Beat; zeros(b,1)];
end