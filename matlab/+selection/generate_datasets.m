function Result = generate_datasets(leadnames, basedir, methodname, ...
    recordmap, discard, countmap, ratiomap, feature)

for i = 1:length(leadnames)
    lead = leadnames{i};
    if countmap(lead) == 0
        continue;
    end
    
    if feature == 0
        disp(['Processing lead ' lead ' for ST segment diagnosis...']);
    else
        disp(['Processing lead ' lead ' for T wave diagnosis...']);
    end
    
    [inputs,targets] = selection.generate_lead_dataset(...
        basedir, lead, methodname, recordmap(lead), ...
        discard, countmap(lead), ratiomap(lead), feature);
    
    Result.(lead) = struct('inputs',inputs,'targets',targets);
end