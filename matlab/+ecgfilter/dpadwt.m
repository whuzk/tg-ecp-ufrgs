function [A,D] = dpadwt(x,J,H,G)

A = cell(J,1);
D = cell(J,1);
for j = 1:J
    a = wconv1(x(:)', H);
    d = wconv1(x(:)', G);
    A{j} = a(2:2:end);
    D{j} = d(2:2:end);
    x = A{j};
end
