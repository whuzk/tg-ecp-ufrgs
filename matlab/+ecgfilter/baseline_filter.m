function Result = baseline_filter(Signal,Fs)
import ecgmath.*

M = round(0.0225*Fs+3.2158);
L1 = 2*M+1;
L2 = 4*M+1;

oc = erosion(dilation(dilation(erosion(Signal,L1),L1),L2),L2);
co = dilation(erosion(erosion(dilation(Signal,L1),L1),L2),L2);
Result = Signal - (oc+co)/2;


function Result = erosion(Signal, M)
N = length(Signal);
Result = zeros(N,1);
delay = (M-1)/2;
pad = zeros(delay,1);
Signal = [pad; Signal; pad];
k = 0:M-1;
for i = M:N+M-1
    Result(i-M+1) = min(Signal(i-k));
end

function Result = dilation(Signal, M)
N = length(Signal);
Result = zeros(N,1);
delay = (M-1)/2;
pad = zeros(delay,1);
Signal = [pad; Signal; pad];
k = 0:M-1;
for i = M:N+M-1
    Result(i-M+1) = max(Signal(i-k));
end