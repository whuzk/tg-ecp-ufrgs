function Result = create_beatset(Database, LeadName)
% Criaçao do conjunto de batidas
import ecgutilities.*

Records = fieldnames(Database);
for i = 1%1:numel(Records)
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
        [Beats,R,RR,Template] = ecgfilter.preprocess(Data(:,j), Fs);
        
        % seleciona as batidas corretamente identificadas
        [index1,index2] = ecgutilities.find_close_beats(Bp, R, Fs);
        
        % agrega anotaçoes de diagnostico
        id = num2str(j-1);
        STseg = extract_diagnosis(ECG.Annotations, Bp(index1), 's', 'ST', id);
        Twave = extract_diagnosis(ECG.Annotations, Bp(index1), 'T', 'T', id);
        
        % armazena o resultado na variavel
        Result.(Records{i}).Fs = Fs;
        Result.(Records{i}).R = R(index2);
        Result.(Records{i}).RR = RR(index2);
        Result.(Records{i}).Beats = Beats(:,index2);
        Result.(Records{i}).Template = Template;
        Result.(Records{i}).STseg = STseg;
        Result.(Records{i}).Twave = Twave;
    end
end
