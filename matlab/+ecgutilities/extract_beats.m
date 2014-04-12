function [Result] = extract_beats(Signal, Rpeaks, RR, FrameSize)

N = length(Signal);
m = length(Rpeaks);

half = floor(min(RR,FrameSize)/2);
Result = zeros(FrameSize,m);
for i = 1:m
    left = Rpeaks(i)-half(i);
    right = Rpeaks(i)+half(i);
    if left < 1
        d = 1-left;
        left = left+d;
        right = right-d;
    elseif right > N
        d = right-N;
        left = left+d;
        right = right-d;
    end
    Beat = Signal(left:right);
    Result(:,i) = frame_beat(Beat,FrameSize);
end

function Result = frame_beat(Signal, L)
D = L-length(Signal);
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    Result = Signal(1+a:end-b);
else
    Result = [zeros(a,1); Signal; zeros(b,1)];
end