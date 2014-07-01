function tables = generate_ind_stats(basedir)
global methods

tables = cell(8,1);
count = 1;
for i = 0:1
    disp(['Training networks for feature ' num2str(i)]);
    for j = 0:1
        disp(['Training networks for channel ' num2str(j)]);
        for k = 1:length(methods)
            if i > 0 && k < 3
                continue;
            end
            disp(['Training networks for method ' methods{k}]);
            tables{count} = get_stat_table(basedir, i, j, methods{k});
            count = count + 1;
        end
    end
end

function Result = get_stat_table(basedir, i, j, method)
global measures
[stats,RecordName] = get_statistics(basedir, i, j, method, 3);
Result = [table(RecordName) array2table(stats, 'VariableNames', measures)];

function [stats,names] = get_statistics(basedir, feature, channel, ...
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
    stats(i,:) = calc_statistics(dataset, feature, max_attempts);
end

function Result = calc_statistics(dataset, feature, max_attempts)
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