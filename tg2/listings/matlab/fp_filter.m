function [sigD,sigM] = fp_filter(x, Fs, delay)
% filtro passa-baixas com frequencia de corte ~11 Hz
[bl,al,~,dl] = intfdesign.lowpass('N,F3db',2,11,Fs);
% filtro de media-movel com largura de ~50 ms
[bi,ai,~,di] = intfdesign.maverage('Width',0.05,Fs);
% filtro 'diferenciador' de 2 pontos
[bd,ad,~,dd] = intfdesign.derivative('N,M',1,0);
% filtro derivativo morfologico de escala s
s = round(0.06*Fs); B = ones(2*s+1,1);
% aplica os filtros
sigL = filter(bl, al, x);
sigI = filter(bi, ai, sigL);
sigD = filter(bd, ad, sigI);
sigM = (imerode(sigL,B)+imdilate(sigL,B)-2*sigL);
% atraso artificial dos sinais
d = ceil(delay - (dl + di + dd));
sigD = [zeros(d,1); sigD(1:end-d)];
d = ceil(delay - dl);
sigM = [zeros(d,1); sigM(1:end-d)];