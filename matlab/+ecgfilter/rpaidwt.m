function a = rpaidwt(a,D,H,G)

A{J} = a;
A = ridwt(0,J-1,A,D,H,G);
for i = 1:N/2-1
    A = ridwt(i,1,A,D,H,G);
end
a = A{1};

function [a,A = ridwt(i,j,A,D,H,G)
if (i == 0)
    if (j == 0)
        A{0+1}(0+1) = A{0+1+1}(0+1)*H(1) + D{j+1}(1)*G(1);
        A{0+1}(1+1) = ?;
    else
        A{j+1}(0+1) = ?;
        A = ridwt(0,j-1,A,D,H,G);
    end
elseif mod(i,2) ~= 0
    A{j+1}((i+1)/2) = ?;
else
    A = ridwt(i/2,j+1,A,D,H,G);
end