function [stats,names] = generate_ind_stats(basedir, feature, channel, ...
    methodname, max_attempts)

files = dir([basedir '*.mat']);
stats = zeros(length(files),6);
names = cell(length(files),1);
for i = 1:length(files)
    file = files(i);
    [~,names{i},~] = fileparts(file.name);
    disp(['Processing ' names{i} '...']);
    filepath = [basedir file.name];
    info = load(filepath);
    fields = fieldnames(info);
    dataset = info.(fields{channel+1}).Datasets.(methodname);
    stats(i,:) = calc_stat(dataset, feature, max_attempts);
end

function Result = calc_stat(dataset, feature, max_attempts)
x = dataset(:,1:end-2)';
t = 2*(dataset(:,end-1+feature)'~=0)-1;
attempts = 0;
minperf = Inf;
while attempts < max_attempts
    layers = randi([5 15],1);
    [net,tr] = nnetwork.train_network(x, t, ...
        layers, 'trainlm', 'tansig', 'mapstd');
    if tr.best_perf < minperf
        minperf = tr.best_perf;
        Result = utils.compute_statistics(t > 0, net(x) > 0);
    end
    attempts = attempts + 1;
end