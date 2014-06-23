function train_gopalak_networks
import utils.*;
import nnetwork.*;
global leadnames

basedir = 'C:\physiobank\database\';
loadpath = [basedir 'gopalak_datasets.mat'];
savepath = [basedir 'gopalak_networks.mat'];
datasets = load(loadpath);
trainfcn = 'trainscg';
transfcn = 'tansig';

for i = 1:length(leadnames)
    lead = leadnames{i};
    if ~isfield(datasets, lead) || isempty(datasets.(lead))
        continue;
    end
    disp(['Generating networks for lead ' lead '...']);
    
    x = datasets.(lead)(:,1:end-2)';
    t1 = 2*(datasets.(lead)(:,end-1)' == 0)-1;
    t2 = 2*(datasets.(lead)(:,end)' == 0)-1;
    layer = randi([10 20],1);
    [net1,tr1] = train_network(x, t1, layer, trainfcn, transfcn);
    compute_statistics(t1 < 0, net1(x) < 0)
    [net2,tr2] = train_network(x, t2, layer, trainfcn, transfcn);
    compute_statistics(t2 < 0, net2(x) < 0)
    networks.(lead).nets = {net1 net2};
    networks.(lead).recs = [tr1 tr2];
end
save(savepath, '-struct', 'networks');