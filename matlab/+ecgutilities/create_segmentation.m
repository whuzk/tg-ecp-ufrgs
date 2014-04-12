function Result = create_segmentation(Database, LeadName)
import ecgutilities.*

Records = fieldnames(Database);
for i = 1:numel(Records)
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
        %continue;
        
        % processa o sinal da derivaçao
        Rpeaks = ecgfeatures.detect_qrs(Data(:,j),Fs);
        [A,B,R] = ecgutilities.merge_rpeaks(Bp, Rpeaks, Fs);
        Result.(Records{i}).A = A;
        Result.(Records{i}).B = B;
        Result.(Records{i}).R = R;
        %{
        figure; hold on;
        plot(Data(:,j));
        stem(R, A*0.2, 'g');
        stem(R, B*0.1, 'r');
        %}
        ecgmath.compute_statistics(A,B)
    end
end