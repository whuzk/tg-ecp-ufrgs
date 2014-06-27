function [ST,T] = get_ratios(count,fields)

N = length(fields);
ST = zeros(N,3);
T = zeros(N,3);
for i = 1:N
    lead = fields{i};
    total = sum(sum(count.(lead)));
    ST(i,:) = sum(count.(lead),2)/total;
    T(i,:) = sum(count.(lead),1)/total;
end