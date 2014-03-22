function [A, B] = find_close_beats(Rpeaks, Predicted)
%   Obtem os incides das batidas que foram corretamente idfentificadas.
%
% Entradas:
%   Rpeaks    - localizaçao correta dos picos da onda R
%   Predicted - localizaçao predita dos picos da onda R
%
% Saída:
%   incides das batidas identificadas
%
n = length(Rpeaks);
m = length(Predicted);
A = false(n,1);
B = false(m,1);

i = 1;
j = 2;
while (i <= m)
    while (j <= n) && (Predicted(i) - Rpeaks(j) > 15)
        j = j + 1;
    end
    if (j > n)
        break
    elseif abs(Predicted(i) - Rpeaks(j)) < 15
        A(j) = true;
        B(i) = true;
        j = j + 1; 
    end
    i = i + 1;
end
