function [a,b,c] = merge_qrs(Rref, Rtest, Fs)

VDI = round(0.15*Fs);

n = length(Rref);
m = length(Rtest);
a = false(n+m,1);
b = false(n+m,1);
c = zeros(n+m,1);

total = max(n,m);
Rref = [Rref(:); inf(1+total-n,1)];
Rtest = [Rtest(:); inf(1+total-m,1)];

i = 1;
j = 1;
k = 1;
while (i <= total) && (j <= total)
    % get differences
    d1 = Rref(i) - Rtest(j);
    d2 = Rref(i+1) - Rtest(j+1);
    
    % check precedence
    if d1 > 0
        % test annotation is earliest
        d3 = Rref(i) - Rtest(j+1);
        
        if (d1 <= VDI && (d1 < abs(d3) || abs(d2) < abs(d3)))
            % (1) current test annotation is within the validation interval
            % and is a better match than the next test annotation: pair it
            a(k) = true;
            b(k) = true;
            c(k) = Rref(i);
            i = i + 1;
            j = j + 1;
        else
            % (2) there is no match to the current test annotation. hence,
            % keep it by itself and go to the next test annotation.
            b(k) = true;
            c(k) = Rtest(j);
            j = j + 1;
        end
    else
        % reference annotation is earliest
        d4 = Rref(i+1) - Rtest(j);
        
        if (-d1 <= VDI && (-d1 < abs(d4) || abs(d2) < abs(d4)))
            % (3) current ref. annotation is within the validation interval
            % and is a better match than the next ref. annotation: pair it
            a(k) = true;
            b(k) = true;
            c(k) = Rref(i);
            j = j + 1;
            i = i + 1;
        else
            % (4) There is no match to the current ref. annotation. hence,
            % keep it by itself and go to the next ref. annotation.
            a(k) = true;
            c(k) = Rref(i);
            i = i + 1;
        end
    end
    
    % advance overall result
    k = k + 1;
end

% trim the output vectors
a(k:end) = [];
b(k:end) = [];
c(k:end) = [];