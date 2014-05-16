function J = mohebbi_jpoints(Beats, defI, defJ, Fs, gain)

[bi,ai,gi,di] = intfdesign.maverage('Width',0.05,Fs);
[bd,ad,~,dd] = intfdesign.derivative('N,M',1,0);
delay = ceil(di + dd);

center = floor(size(Beats,1)/2)+1;
L1 = round(0.02*Fs);
L2 = round(0.12*Fs);
thr = gi*gain*1.25/Fs;
M = size(Beats,2);
J = zeros(M,1);
for i = 1:M
    sigI = filter(bi, ai, Beats(:,i));
    sigD = filter(bd, ad, sigI);
    istart = center + L1 + delay;
    iend = center + L2 + delay;
    J(i) = get_point(abs(sigD), istart, iend, L1, thr, defJ(i)) - delay;
end


function pos = get_point(data, istart, iend, L, thr, default)
M = floor(L/2);
pos = default;
for i = istart+L-1:iend
    if all(data(i-L+1:i) <= thr)
        pos = i - M;
        break;
    end
end