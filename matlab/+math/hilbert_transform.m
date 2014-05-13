function Result = hilbert_transform(x)

n = length(x);
h = zeros(n,1);
if mod(n,2) == 0
  % even
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
else
  % odd
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end

Result = ifft(fft(x).*h);