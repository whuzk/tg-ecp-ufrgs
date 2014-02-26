function Result = create_beatset(Database, LeadName)
% Cria�ao do conjunto de batidas
%
import ecgutilities.*;

Records = fieldnames(Database);
for i = 4:4%numel(Records)
    ECG = Database.(Records{i});
    
    [Fs,Bp,Sc,Leads,Data] = ecg_interpret(ECG);
    for j = 1:Sc
        % seleciona a deriva�ao correta
        if ~strcmp(Leads{j}, LeadName)
            disp(['ignoring ' Records{i} '.' Leads{j}]);
            continue;
        else
            disp(['processing ' Records{i} '.' Leads{j}]);
        end
        
        % processa o sinal da deriva�ao
        [Beats,Rpeaks,RR,Template] = ecg_preprocess(Data(:,j), Fs);
        
        % seleciona as batidas corretamente identificadas
        [index1,index2] = ecg_find_close_beats(Bp, Rpeaks-5-2);
        
        % agrega anota�oes de diagnostico
        id = num2str(j-1);
        STseg = ecg_get_diagnosis(ECG.Annotations, Bp(index1), 's', 'ST', id);
        Twave = ecg_get_diagnosis(ECG.Annotations, Bp(index1), 'T', 'T', id);
        
        % armazena o resultado na variavel
        Result.(Records{i}).Fs = Fs;
        Result.(Records{i}).RR = RR(index2);
        Result.(Records{i}).Beats = Beats(:,index2);
        Result.(Records{i}).Template = Template;
        Result.(Records{i}).STseg = STseg;
        Result.(Records{i}).Twave = Twave;
    end
end
