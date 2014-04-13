function Y = wavelet_baseline(X,Fs,wname)

% initializations
N = length(X);
J = round(log2(Fs));

% load wavelet filters
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);

% wavelet decomposition
[A,D] = ecgwavelet.dpadwt2(X,J,Lo_D,Hi_D);

% remove high-frequency noise
for j = 1:J
    D{j} = zeros(size(D{j}));
end

% wavelet reconstruction
Y = ecgwavelet.dpaidwt2(A{end},D,N,Lo_R,Hi_R);