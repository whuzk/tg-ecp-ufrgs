function test_segmentation(Database, LeadName)
import ecgutilities.*

Records = fieldnames(Database);
M = numel(Records);
A = cell(M,1);
B = cell(M,1);
count = 0;
for i = 1:M
    ECG = Database.(Records{i});
    
    [Fs,Bp,Sc,Leads,Data] = ecgutilities.interpret(ECG);
    for j = 1:Sc
        % seleciona a derivaçao correta
        if ~strcmp(Leads{j}, LeadName)
            disp(['ignoring ' Records{i} '.' Leads{j}]);
            continue;
        else
            disp(['processing ' Records{i} '.' Leads{j}]);
        end
        
        % processa o sinal da derivaçao
        Rpeaks = ecgfilter.detect_qrs(Data(:,j),Fs);
        [A{i},B{i}] = ecgutilities.merge_rpeaks(Bp, Rpeaks);
        count = count + length(A{i});
        
        %{
        figure; hold on;
        plot(Data(:,j));
        stem(R, A{i}*0.2, 'g');
        stem(R, B{i}*0.1, 'r');
        %}
    end
end

%% agrupa os resultados
A_all = false(count,1);
B_all = false(count,1);
iStart = 1;
iEnd = 0;
for i = 1:M
    iEnd = iEnd + length(A{i});
    A_all(iStart:iEnd) = A{i};
    B_all(iStart:iEnd) = B{i};
    iStart = iEnd + 1;
end
ecgmath.compute_statistics(A_all,B_all)
