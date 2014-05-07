function [sigD,sigM,gainD,delay] = yansun_filter(x, Fs)
% Funçao para transformar o sinal de ECG num sinal que possa ser utilizado
% pelo algoritmo de detecçao de batimentos

% filtro passa-baixas com frequencia de corte ~11 Hz
[bl,al,gl,dl] = intfdesign.lowpass('N,F3db',2,11,Fs);

% filtro de media-movel com largura de ~20 ms
[bi,ai,gi,di] = intfdesign.maverage('Width',0.02,Fs);

% filtro derivativo de 2 pontos
[bd,ad,gd,dd] = intfdesign.derivative('N,M',1,0);

% filtro derivativo morfologico de escala s
s = round(0.06*Fs);
B = ones(2*s+1,1);

% aplica os filtros
sigL = filter(bl, al, x);   %./(2^nextpow2(gi));%./gl;
sigI = filter(bi, ai, sigL);%./(2^nextpow2(gi));%./gl;
sigD = filter(bd, ad, sigI);%./(2^nextpow2(gd));%./gd;
sigM = (imerode(sigL,B)+imdilate(sigL,B)-2*sigL);%./s;

% atrasa os sinais
d = s - ceil(di + dd);
sigD = [zeros(d,1); sigD(1:end-d)];
sigM = [zeros(s,1); sigM(1:end-s)];

% calcula o atraso total e o ganho do sinal derivativo
gainD = gl * gi;
delay = dl + s;