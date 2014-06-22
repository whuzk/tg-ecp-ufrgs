function Result = generate_datasets(basedir, methodname, ...
    recordmap, discard, countmap, ratio)
global leadnames

for i = 1:length(leadnames)
    lead = leadnames{i};
    if countmap(lead) == 0
        continue;
    end
    disp(['Processing lead ' lead '...']);
    Result.(lead) = nnetwork.generate_lead_dataset(basedir, lead, ...
        methodname, recordmap(lead), discard, countmap(lead), ratio);
end