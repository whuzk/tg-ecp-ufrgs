function Result = generate_datasets(leadnames, basedir, methodname, ...
    recordmap, discard, countmap, stratiomap, tratiomap)

for i = 1:length(leadnames)
    lead = leadnames{i};
    
    disp(['Processing lead ' lead ' for ST segment diagnosis...']);
    [inputs,targets] = nnetwork.generate_lead_dataset(...
        basedir, lead, methodname, recordmap(lead), ...
        discard, countmap(lead), stratiomap(lead), false);
    Result.(lead).ST = struct('inputs',inputs,'targets',targets);
    
    disp(['Processing lead ' lead ' for T wave diagnosis...']);
    [inputs,targets] = nnetwork.generate_lead_dataset(...
        basedir, lead, methodname, recordmap(lead), ...
        discard, countmap(lead), tratiomap(lead), true);
    Result.(lead).T = struct('inputs',inputs,'targets',targets);
end