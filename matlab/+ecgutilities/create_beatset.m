function Result = create_beatset(Database)
% Criaçao do conjunto de batidas
import ecgutilities.*

Records = fieldnames(Database);
for i = 1%1:numel(Records)
    ECG = Database.(Records{i});
    [Fs,Bp,Sc,Leads] = ecgutilities.interpret(ECG);
    for j = 1%1:Sc
        disp(['processing ' Records{i} '.' Leads{j}.name]);
        Info = extract_beat_info(Leads{j}.data,Fs,Bp,num2str(j-1),ECG.Annotations);
        Result.(Records{i}).(Leads{j}.name) = Info;
    end
end
