function Result = count_classes(basedir)

files = dir([basedir '*.mat']);
Result = struct;

for i = 1:length(files)
    file = files(i);
    if isempty(file.name) || file.isdir
        continue;
    end
    disp(['Processing ' file.name '...']);
    
    record = load([basedir file.name]);
    for j = 1:record.SignalCount
        lead = record.Info(j).Description;
        atr = record.Annotations.atr;
        idx = ismember(atr.Type, 'NLRBAaJSVFejnE/fQ?');
        qrs = record.Annotations.atr.Sample(idx);
        ST = utils.extract_diagnosis(qrs, atr, 's', 'ST', j-1);
        T = utils.extract_diagnosis(qrs, atr, 'T', 'T', j-1);
        
        stnormal = ST.diagVal == 0;
        stelevated = ST.diagVal > 0;
        stdepressed = ST.diagVal < 0;
        tnormal = T.diagVal == 0;
        televated = T.diagVal > 0;
        tdepressed = T.diagVal < 0;
        
        C = zeros(3,3);
        C(1,1) = C(1,1) + length(find(stnormal & tnormal));
        C(2,1) = C(2,1) + length(find(stelevated & tnormal));
        C(3,1) = C(3,1) + length(find(stdepressed & tnormal));
        C(1,2) = C(1,2) + length(find(stnormal & televated));
        C(2,2) = C(2,2) + length(find(stelevated & televated));
        C(3,2) = C(3,2) + length(find(stdepressed & televated));
        C(1,3) = C(1,3) + length(find(stnormal & tdepressed));
        C(2,3) = C(2,3) + length(find(stelevated & tdepressed));
        C(3,3) = C(3,3) + length(find(stdepressed & tdepressed));
        
        if ~isfield(Result, lead)
            Result.(lead) = C;
        else
            Result.(lead) = Result.(lead) + C;
        end
    end
end