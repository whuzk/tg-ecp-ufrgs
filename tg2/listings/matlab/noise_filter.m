function sigN = noise_filter(x, Fs, Fm, delay)
% filtro rejeita-faixa com frequencia central em Fm
[bs,as] = cheby1(2,0.1,[Fm-1 Fm+1]*2/Fs,'stop');
ds = mean(grpdelay(bs,as,[5 15]*2/Fs));
% filtro passa-baixas com frequencia de corte ~40Hz
[bl,al] = butter(4,40*2/Fs);
dl = mean(grpdelay(bl,al,[5 15]*2/Fs));
% aplica os filtros
sigS = filter(bs, as, x);
sigL = filter(bl, al, sigS);
% atraso artificial so sinal
d = ceil(delay - ds - dl);
sigL = [zeros(d,1); sigL(1:end-d)];