function Result = ecg_get_diagnosis(Annot, Rpeaks, Type, keyWord, signalID)

m = length(Rpeaks);
l = length(keyWord);
B = str2double(Annot.Sample);
index = find(strcmp(Annot.Type, Type));

Result.elevation = false(m,1);
Result.depression = false(m,1);
Result.elevPeakValue = zeros(m,1);
Result.depPeakValue = zeros(m,1);

k = 1;
while k <= length(index)
    aux = Annot.Aux{index(k)};
    if (~isempty(aux) && aux(1) == '(' && ...
        strcmp(aux(1+(1:l)), keyWord) && aux(2+l) == signalID)
        % inicio de episodio
        beginIndex = index(k);
        ch = aux(3+l);
        
        % procura a marcaçao de pico
        k = k + 1;
        while k <= length(index)
            aux = Annot.Aux{index(k)};
            if (~isempty(aux) && upper(aux(1)) == 'A' && ...
                strcmp(aux(1+(1:l)), keyWord) && aux(2+l) == signalID)
                % fim de episodio
                peakDev = str2double(aux(4+l:end))/1000;
                break;
            else
                k = k + 1;
            end
        end
        
        % procura a marcaçao de termino
        k = k + 1;
        while k <= length(index)
            aux = Annot.Aux{index(k)};
            if (~isempty(aux) && aux(end) == ')' && ...
                strcmp(aux(1:l), keyWord) && aux(1+l) == signalID)
                % fim de episodio
                endIndex = index(k);
                break;
            else
                k = k + 1;
            end
        end
        
        % verifica se é elevaçao ou depressao
        if k <= length(index)
            if ch == '+'
                j = B(beginIndex) <= Rpeaks & Rpeaks <= B(endIndex);
                Result.elevation(j) = true;
                Result.elevPeakValue(j) = peakDev;
            elseif ch == '-'
                j = B(beginIndex) <= Rpeaks & Rpeaks <= B(endIndex);
                Result.depression(j) = true;
                Result.depPeakValue(j) = peakDev;
            end
        end
    end
    k = k + 1;
end
