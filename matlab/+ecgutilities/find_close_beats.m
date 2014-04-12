function [A, B] = find_close_beats(R1, R2, Fs)
%   Obtem os incides das batidas que foram corretamente idfentificadas.
%
% Entradas:
%   R1 - localizaçao correta dos picos da onda R
%   R2 - localizaçao predita dos picos da onda R
%
% Saída:
%   incides das batidas identificadas
%
VDI1 = round(0.10*Fs);
VDI2 = round(0.20*Fs);

n = length(R1);
m = length(R2);
A = false(n,1);
B = false(m,1);

i = 1;
j = 1;
while (i <= n) && (j <= m)
    while (i <= n) && (j <= m) && (R2(j) <= R1(i))
        if R1(i) - R2(j) <= VDI1
            A(i) = true;
            B(j) = true;
            i = i + 1;
        end
        j = j + 1;
    end
    while (i <= n) && (j <= m) && (R1(i) <= R2(j))
        if R2(j) - R1(i) <= VDI2
            A(i) = true;
            B(j) = true;
            j = j + 1;
        end
        i = i + 1;
    end
end