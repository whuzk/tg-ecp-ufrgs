function [sigI,delay] = qrs_filter(x, Fs)
% filtro passa-baixas com frequencia de corte ~11 Hz
[bl,al,gl,dl] = intfdesign.lowpass('N,F3db',2,11,Fs);
% filtro passa-altas com frequencia de corte ~5 Hz
[bh,ah,gh,dh] = intfdesign.highpass('N,F3db',1,5,Fs);
% filtro 'diferenciador' de 5 pontos
[bd,ad,gd,dd] = intfdesign.derivative('N,M',1,3);
% filtro de media-movel com largura de ~150 ms
[bi,ai,gi,di] = intfdesign.maverage('Width',0.15,Fs);
% aplica os filtros
sigL = filter(bl, al, x)    ./gl;
sigF = filter(bh, ah, sigL) ./gh;
sigD = filter(bd, ad, sigF) ./gd;
sigS = min(sigD.^2./2^3, 2^15-1);
sigI = filter(bi, ai, sigS) ./gi;
% calcula o atraso total
delay = dl + dh + dd + di;