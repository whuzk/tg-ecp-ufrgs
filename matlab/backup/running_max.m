function [Result,mem,cmp] = running_max(Signal,M,ai)
N = length(Signal);
Result = zeros(N,1);
mem = zeros(N,1);
cmp = zeros(N,1);
pos = zeros(1,M);
val = zeros(1,M);
first = 1;
len = 1;
pos(1) = 0;
val(1) = ai;
for i = 1:N
    % get the current sample
    a = Signal(i);
    
    % search for the first element greater than the current sample
    j = len;
    count = 1;
    while j > 0 && a >= val(mod(first+j-2,M)+1)
        j = j - 1;
        count = count + 1;
    end
    
    % put the sample next to element found and adjust the length
    idx = mod(first+j-1,M)+1;
    val(idx) = a;
    pos(idx) = i;
    len = j + 1;
    
    % check if the first in line has gone out of the windows length
    if len > M || pos(first) <= i-M
        len = len - 1;
        first = mod(first,M)+1;
    end
    
    % store the results
    Result(i) = val(first);
    mem(i) = len;
    cmp(i) = count;
end