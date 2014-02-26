function [x,y,g] = ecg_baseline_points(Beat, Fs, R)

slope = 2.5/Fs;
g = gradient(Beat);
%GG = smooth(g,5);
Gabs = abs(g);

N = length(Beat);
l = ceil(0.01*Fs);
left = max(1,R-fix(0.08*Fs));
right = min(N,R-fix(0.04*Fs));

i1 = find(Gabs <= slope, 1, 'first')+l;
if isempty(i1)
    i1 = 1+l;
end

i2 = left+find(Gabs(left:right) <= slope, 1, 'first')-1+l;
if isempty(i2)
    i2 = left+l;
end

i3 = find(Gabs <= slope, 1, 'last')-l;
if isempty(i3)
    i3 = length(Beat)-l;
end

i_aux = i1+fix((i2-i1)*[1/5 2/5 3/5]);

x = unique([i1 i_aux i2 i3]);
y = zeros(1,length(x)+2);
y(1) = mean(g(x(1)-l:x(1)+l));
for j = 1:length(x)
    y(j+1) = mean(Beat(x(j)-l:x(j)+l));
end
y(end) = 0;

%figure, plot(g), grid on;
%figure, plot(GG), grid on;
%figure, plot(Beat), grid on;
