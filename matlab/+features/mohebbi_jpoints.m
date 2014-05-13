function J = mohebbi_jpoints(Beats, R, default, Fs, gain)

[bi,ai,gi,di] = intfdesign.maverage('Width',0.05,Fs);
[bd,ad,~,dd] = intfdesign.derivative('N,M',1,0);
delay = ceil(di + dd);

L1 = round(0.02*Fs);
L2 = round(0.12*Fs);
thr = gain*1.25/Fs;
M = size(Beats,2);
J = zeros(M,1);
for i = 1:M
    sigI = filter(bi, ai, Beats(:,i)) ./ gi;
    sigD = filter(bd, ad, sigI);
    istart = R + L1 + delay;
    iend = R + L2 + delay;
    J(i) = get_point(abs(sigD), istart, iend, L1, thr, default(i)) - delay;
end


function pos = get_point(data, istart, iend, L, thr, default)
M = floor(L/2);
pos = default;
for i = istart+M:iend-M
    if all(data(i-M:i+M) < thr)
        pos = i;
        break;
    end
end