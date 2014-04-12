function Result = suppress_baseline(Beat,l)

n = length(Beat);
y0 = mean(Beat(1:l));
y1 = mean(Beat(end-l+1:end));
P = polyfit([1 n],[y0 y1],1);
Result = Beat-polyval(P,1:n)';