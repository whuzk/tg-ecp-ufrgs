function Y = wavfilter(X,J,wname)
import ecgfilter.*

N = length(X);

% load wavelet filters
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);

% wavelet decomposition
[A,D] = dpadwt2(X,J,Lo_D,Hi_D);

% remove low-frequency noise
a0 = zeros(size(A{end}));

% remove high-frequency noise
D = denoise(D,J);

% wavelet reconstruction
Y = dpaidwt2(a0,D,N,Lo_R,Hi_R);


function D = denoise(D,J)
thr = 4*std(D{1});
for j = 1:J
    D{j} = D{j}.*(abs(D{j}) >= thr);
end