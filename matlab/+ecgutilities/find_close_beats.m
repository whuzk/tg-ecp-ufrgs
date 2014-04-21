function [A, B] = find_close_beats(Rref, Rtest, Fs)
%   Obtem os incides das batidas que foram corretamente idfentificadas.
%
% Entradas:
%   R1 - localizaçao correta dos picos da onda R
%   R2 - localizaçao predita dos picos da onda R
%
% Saída:
%   incides das batidas identificadas
%
VDI = round(0.15*Fs);

n = length(Rref);
m = length(Rtest);
A = false(n,1);
B = false(m,1);

Rref = [Rref(:); Inf];
Rtest = [Rtest(:); Inf];

i = 1;
j = 1;
while (i <= n) && (j <= m)
    % get differences
    d1 = Rref(i) - Rtest(j);
    d2 = Rref(i+1) - Rtest(j+1);
    
    % check precedence
    if d1 > 0
        % test annotation is earliest
        d3 = Rref(i) - Rtest(j+1);
        
        if (d1 <= VDI && (d1 < abs(d3) || abs(d2) < abs(d3)))
            % (1) current test annotation is within the validation interval
            % and is a better match than the next test annotation: pair it
            A(i) = true;
            B(j) = true;
            i = i + 1;
            j = j + 1;
        else
            % (2) there is no match to the current test annotation. hence,
            % do not do anything and go to the next test annotation.
            j = j + 1;
        end
    else
        % reference annotation is earliest
        d4 = Rref(i+1) - Rtest(j);
        
        if (-d1 <= VDI && (-d1 < abs(d4) || abs(d2) < abs(d4)))
            % (3) current ref. annotation is within the validation interval
            % and is a better match than the next ref. annotation: pair it
            A(i) = true;
            B(j) = true;
            j = j + 1;
            i = i + 1;
        else
            % (4) There is no match to the current ref. annotation. hence,
            % do not do anything and go to the next ref. annotation.
            i = i + 1;
        end
    end
end