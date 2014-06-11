function y = running_minmax(x,M,ismax)
N = length(x);
y = zeros(N,1);
pos = zeros(1,M);
val = zeros(1,M);
first = 1;
count = 1;
for i = 1:N
    % get the current sample
    a = x(i);
    
    % search for the first element greater than the current sample
    j = count;
    if ismax
        while j > 0 && a >= val(mod(first+j-2,M)+1)
            j = j - 1;
        end
    else
        while j > 0 && a <= val(mod(first+j-2,M)+1)
            j = j - 1;
        end
    end
    
    % put the sample next to element found and adjust the length
    idx = mod(first+j-1,M)+1;
    val(idx) = a;
    pos(idx) = i;
    count = j + 1;
    
    % check if the first in line has gone out of the windows length
    if count > M || pos(first) <= i-M
        count = count - 1;
        first = mod(first,M)+1;
    end
    
    % store the results
    y(i) = val(first);
end