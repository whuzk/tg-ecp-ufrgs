function Result = frame_beat(Beat, R, L)

Result = zeros(L,1);
half = floor(L/2);
N = length(Beat);

trim = max(0,R-1-half);
pad = max(0,half-R+1);
Result(1:half+1) = [zeros(pad,1); Beat(trim+1:R)];

trim = max(0,N-R-half);
pad = max(0,half-N+R);
Result(half+2:end) = [Beat(R+1:end-trim); zeros(pad,1)];