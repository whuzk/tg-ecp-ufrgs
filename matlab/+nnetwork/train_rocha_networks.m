function train_rocha_networks
import utils.*;
import nnetwork.*;
global leadnames

basedir = 'C:\physiobank\database\';
loadpath = [basedir 'rocha_datasets.mat'];
savepath = [basedir 'rocha_networks.mat'];
datasets = load(loadpath);
layersmap = containers.Map(leadnames, ...
    {[] [] [5 3] [3 3] [] [5 3] [] [] [] [5 3] [5 3] [5 3] [6 4] [7 3] []});
trainfcn = 'trainlm';
transfcn = 'tansig';

for i = 1:length(leadnames)
    lead = leadnames{i};
    if ~isfield(datasets, lead)
        continue;
    end
    disp(['Generating networks for lead ' lead '...']);
    
    x = datasets.(lead)(:,1:end-2)';
    t1 = 2*(datasets.(lead)(:,end-1)' > 0)-1;
    t2 = 2*(datasets.(lead)(:,end-1)' < 0)-1;
    [net1,tr1] = train_network(x, t1, layersmap(lead), trainfcn, transfcn);
    [net2,tr2] = train_network(x, t2, layersmap(lead), trainfcn, transfcn);
    compute_statistics(t1 > 0, net1(x) > 0)
    compute_statistics(t2 > 0, net2(x) > 0)
    networks.(lead).nets = [net1 net2];
    networks.(lead).recs = [tr1 tr2];
end
save(savepath, '-struct', 'networks');