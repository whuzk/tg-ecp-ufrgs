function Result = running_std(Signal,M)
N = length(Signal);
Result = zeros(N,1);
mean = 0;
n = 0;
M2 = 0;
for i = 1:M
    [M2,n,mean] = add_variable(M2,n,mean,Signal(i));
    Result(i) = sqrt(M2/(n-1));
end
for i = M+1:N
    [M2,mean] = update_variable(M2,M,mean,Signal(i-M),Signal(i));
    Result(i) = sqrt(M2/(M-1));
end

function [M2,n,mean] = add_variable(M2,n,mean,x)
n = n+1;
delta = x-mean;
mean = mean+delta/n;
M2 = M2+delta*(x-mean);
 
function [M2,n,mean] = remove_variable(M2,n,mean,x)
n = n-1;
delta = x-mean;
mean = mean-delta/n;
M2 = M2-delta*(x-mean);

function [M2,mean] = update_variable(M2,n,mean,oldX,newX)
delta = newX-oldX;
dold = oldX-mean;
mean = mean+delta/n;
dnew = newX-mean;
M2 = M2+delta*(dold+dnew);