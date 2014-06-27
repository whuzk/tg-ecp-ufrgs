function Result = train_networks(leadnames, datasets, trainfcn, transfcn, processfcn)

for i = 1:length(leadnames)
    lead = leadnames{i};
    disp(['Generating networks for lead ' lead '...']);
    [net1,tr1] = inner_train(datasets, lead, 'ST', trainfcn, transfcn, processfcn);
    [net2,tr2] = inner_train(datasets, lead, 'T', trainfcn, transfcn, processfcn);
    Result.nets = {net1 net2};
    Result.recs = {tr1 tr2};
end

function [net,tr] = inner_train(datasets, lead, feature, trainfcn, ...
    transfcn, processfcn)
layers = randi([10 20],1);
x = datasets.(lead).(feature).inputs;
t = datasets.(lead).(feature).targets;
[net,tr] = nnetwork.train_network(x, t, layers, trainfcn, transfcn, processfcn);
[~,i] = max(net(x),[],1);
utils.compute_statistics(~t(1,:), i > 1)