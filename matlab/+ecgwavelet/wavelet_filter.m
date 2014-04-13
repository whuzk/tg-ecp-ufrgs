function Y = wavelet_filter(X,Fs,wname)

% initializations
N = length(X);
J = round(log2(Fs));

% load wavelet filters
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);

% wavelet decomposition
[A,D] = ecgwavelet.dpadwt2(X,J,Lo_D,Hi_D);

% remove low-frequency noise
a0 = zeros(size(A{end}));

% remove high-frequency noise
D = ecgwavelet.wavelet_denoise(D,J,4);

% wavelet reconstruction
Y = ecgwavelet.dpaidwt2(a0,D,N,Lo_R,Hi_R);