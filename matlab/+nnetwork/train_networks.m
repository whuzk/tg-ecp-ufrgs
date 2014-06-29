function Result = train_networks(leadnames, datasets, trainfcn, transfcn, processfcn)

for i = 1:length(leadnames)
    lead = leadnames{i};
    if ~isfield(datasets, lead)
        continue;
    end
    disp(['Generating networks for lead ' lead '...']);
    [net,tr,stat] = inner_train(datasets, lead, trainfcn, transfcn, processfcn);
    Result.(lead).net = net;
    Result.(lead).tr = tr;
    Result.(lead).stat = stat;
end

function [net,tr,stat] = inner_train(datasets, lead, trainfcn, transfcn, processfcn)
layers = randi([10 20],1);
x = datasets.(lead).inputs;
t = datasets.(lead).targets;
[net,tr] = nnetwork.train_network(x, t, layers, trainfcn, transfcn, processfcn);
stat = utils.compute_statistics(t > 0, net(x) > 0);