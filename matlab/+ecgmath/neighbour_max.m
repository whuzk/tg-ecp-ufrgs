function Points = neighbour_max(Signal,Points,Nsize)

N = length(Signal);
for i = 1:length(Points)
    left = max(1,Points(i)-Nsize);
    right = min(N,Points(i)+Nsize);
    [~,x] = max(Signal(left:right));
    Points(i) = left + x - 1;
end