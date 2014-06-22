function train_mohebbi_networks
import utils.*;
import nnetwork.*;
global leadnames

basedir = 'C:\physiobank\database\';
loadpath = [basedir 'mohebbi_datasets.mat'];
savepath = [basedir 'mohebbi_networks.mat'];
datasets = load(loadpath);
trainfcn = 'traingdx';
transfcn = 'logsig';

for i = 1:length(leadnames)
    lead = leadnames{i};
    if ~isfield(datasets, lead) || isempty(datasets.(lead))
        continue;
    end
    disp(['Generating networks for lead ' lead '...']);
    
    x = datasets.(lead)(:,1:end-2)';
    t = datasets.(lead)(:,end-1)' ~= 0 | datasets.(lead)(:,end)' ~= 0;
    [net,tr] = train_network(x, [t; ~t], 20, trainfcn, transfcn);
    y = net(x);
    compute_statistics(t, y(1,:) > y(2,:))
    networks.(lead).nets = net;
    networks.(lead).recs = tr;
end
save(savepath, '-struct', 'networks');