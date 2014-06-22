function train_networks(loadpath, savepath, trainfcn)
import nnetwork.*;

networks = struct;
datasets = load(loadpath);
fields = fieldnames(datasets);
for i = 1:length(fields)
    lead = fields{i};
    if isempty(datasets.(lead))
        continue;
    end
    disp(['Generating networks for lead ' lead '...']);
    layers = [randi([5 15],1) randi([5 15],1)];
    networks.(lead) = train_network(datasets.(lead), layers, trainfcn);
end
save(savepath, '-struct', 'networks');