function Result = extract_diagnosis(Rk, atr, type, keyWord, id)

isqtable = atr(atr.Type == type, {'Sample', 'Comment'});

M = length(Rk);
N = height(isqtable);
L = length(keyWord);

Result.elevation = false(M,1);
Result.depression = false(M,1);
Result.elevPeakValue = zeros(M,1);
Result.depPeakValue = zeros(M,1);

k = 1;
while k <= N
    aux = isqtable.Comment{k};
    if (~isempty(aux) && ...
            aux(1) == '(' && ...
            strcmp(aux(1+(1:L)), keyWord) && ...
            str2double(aux(2+L)) == id)
        % inicio de episodio
        beginIndex = isqtable.Sample(k);
        ch = aux(3+L);
        
        % procura a marca�ao de pico
        k = k + 1;
        while k <= N
            aux = isqtable.Comment{k};
            if (~isempty(aux) && ...
                    upper(aux(1)) == 'A' && ...
                    strcmp(aux(1+(1:L)), keyWord) && ...
                    str2double(aux(2+L)) == id)
                % fim de episodio
                peakDev = str2double(aux(4+L:end));
                break;
            else
                k = k + 1;
            end
        end
        
        % procura a marca�ao de termino
        k = k + 1;
        while k <= N
            aux = isqtable.Comment{k};
            if (~isempty(aux) && ...
                    aux(end) == ')' && ...
                    strcmp(aux(1:L), keyWord) && ...
                    str2double(aux(1+L)) == id)
                % fim de episodio
                endIndex = isqtable.Sample(k);
                break;
            else
                k = k + 1;
            end
        end
        
        % verifica se � eleva�ao ou depressao
        if k <= N
            if ch == '+'
                j = beginIndex <= Rk & Rk <= endIndex;
                Result.elevation(j) = true;
                Result.elevPeakValue(j) = peakDev;
            elseif ch == '-'
                j = beginIndex <= Rk & Rk <= endIndex;
                Result.depression(j) = true;
                Result.depPeakValue(j) = peakDev;
            end
        end
    end
    k = k + 1;
end
