function Result = get_h(b,a,w)

s = exp(1i*w*pi);
Result = abs(horner(b,s) ./ horner(a,s));


function y = horner(p,x)
y = zeros(size(x));
nc = length(p);
if nc > 0
    y(:) = p(1);
end
for i = 2:nc
    y = x .* y + p(i);
end