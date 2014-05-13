function [Beats,F] = extract_beats(data,F,Fs)

N = length(data);
half = floor(Fs*0.6);
FrameSize = 2*half+1;
center = half + 1;

M = 10;%size(F.R,1);
Beats = zeros(FrameSize,M);
for i = 1:M
    Beat = data(max(1,F(i,1)-5):min(N,F(i,7)+5));
    Beat = Beat - polyfit_baseline(Beat,5);
    r = F(i,4) - F(i,1) + 5 + 1;
    Beat = frame_beat(Beat,r,FrameSize);
    F(i,1:3) = max(1,center - (F(i,4) - F(i,1:3)));
    F(i,5:7) = min(FrameSize,center + (F(i,5:7) - F(i,4)));
    F(i,4) = center;
    Beats(:,i) = Beat;
end

function Result = polyfit_baseline(Beat,l)
n = length(Beat);
y0 = mean(Beat(1:l));
y1 = mean(Beat(end-l+1:end));
P = polyfit([0 n-1],[y0 y1],1);
Result = polyval(P,(0:n-1)');

function Result = frame_beat(Beat, r, L)
half = floor(L/2);
N = length(Beat);
trim1 = max(0,r-1-half);
pad1 = max(0,half-r+1);
trim2 = max(0,N-r-half);
pad2 = max(0,half-N+r);
Result = [zeros(pad1,1); Beat(trim1+1:end-trim2); zeros(pad2,1)];