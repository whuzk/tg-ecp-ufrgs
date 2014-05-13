function sigL = noise_filter(x, Fs, delay)

if mod(Fs,50) == 0
    Fm = 50;
else
    Fm = 60;
end

% filtro rejeita-banda com frequencia de corte na freq. da rede eletrica
[bs,as] = cheby1(2,0.1,[Fm-1 Fm+1]*2/Fs,'stop');
ds = mean(grpdelay(bs,as,[5 15]*2/Fs));

% filtro passa-baixas com frequencia de corte ~40Hz
[bl,al] = butter(4,40*2/Fs);
dl = mean(grpdelay(bl,al,[5 15]*2/Fs));

% aplica os filtros
sigS = filter(bs, as, x);
sigL = filter(bl, al, sigS);

% atrasa o sinal
d = ceil(delay - ds - dl);
sigL = [zeros(d,1); sigL(1:end-d)];