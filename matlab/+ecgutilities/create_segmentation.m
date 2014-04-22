function Result = create_segmentation(Database)
import ecgutilities.*

Records = fieldnames(Database);
for i = 1:numel(Records)
    ECG = Database.(Records{i});
    [Fs,Bp,Sc,Leads] = interpret(ECG);
    for j = 1:Sc
        disp(['processing ' Records{i} '.' Leads{j}.name]);
        Info = extract_qrs_info(Leads{j}.data,Fs,Bp);
        Result.(Records{i}).(Leads{j}.name) = Info;
    end
end