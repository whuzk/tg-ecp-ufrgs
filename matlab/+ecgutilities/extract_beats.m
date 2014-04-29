function Result = extract_beats(data,lead,Fs,R,RR,FrameSize)
import ecgutilities.*
import ecgfilter.*

N = length(data);
L = floor(min(RR,FrameSize)/2);
if L(1) > R(1)-1
    L(1) = R(1)-1;
end
if L(end) > N-R(end)
    L(end) = N-R(end);
end

M = length(R);
Result = zeros(FrameSize,M);
for i = 1:1000%M
    Temp = data(R(i)-L(i):R(i)+L(i));
    F = fiducial_marks(Temp,lead,L(i)+1,Fs);
    Temp = suppress_baseline(Temp(F.P(1):F.T(3)),5);
    Result(:,i) = frame_beat(Temp,FrameSize);
end