function [sigI,delay] = tompkins_filter(x, Fs)
% Funçao para transformar o sinal de ECG num sinal que possa ser utilizado
% pelo algoritmo de detecçao de batimentos

% filtro passa-baixas com frequencia de corte ~11 Hz
[bl,al,gl,dl] = intfdesign.lowpass('N,F3db',2,11,Fs);

% filtro passa-altas com frequencia de corte ~5 Hz
[bh,ah,gh,dh] = intfdesign.highpass('N,F3db',1,5,Fs);

% filtro derivativo de 5 pontos
[bd,ad,gd,dd] = intfdesign.derivative('N,M',1,3);

% filtro de media-movel com largura de ~150 ms
[bi,ai,gi,di] = intfdesign.maverage('Width',0.15,Fs);

% aplica os filtros
sigL = filter(bl, al, x)    ./ (2^nextpow2(gl));%./gl;
sigF = filter(bh, ah, sigL) ./ (2^nextpow2(gh));%./gh;
sigD = filter(bd, ad, sigF);%./(2^nextpow2(gd));%./gd;
sigS = min(sigD.^2 ./ 2^3, 2^15-1);
sigI = filter(bi, ai, sigS);%./(2^nextpow2(gi))%./gi;

% calcula o atraso total
delay = dl + dh + dd + di;