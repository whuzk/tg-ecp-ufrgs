import nnetwork.*
global edbleadnames
basedir = 'C:\physiobank\database\';

%%
loadpath = [basedir 'rocha_datasets.mat'];
savepath = [basedir 'rocha_networks.mat'];
netowrks = train_networks(edbleadnames, ...
    load(loadpath), 'trainlm', 'logsig', 'mapstd');
save(savepath, '-struct', 'networks');

%%
loadpath = [basedir 'mohebbi_datasets.mat'];
savepath = [basedir 'mohebbi_networks.mat'];
netowrks = train_networks(edbleadnames, ...
    load(loadpath), 'trainlm', 'logsig', 'mapstd');
save(savepath, '-struct', 'networks');

%%
loadpath = [basedir 'rocha_datasets.mat'];
savepath = [basedir 'rocha_networks.mat'];
netowrks = train_networks(edbleadnames, ...
    load(loadpath), 'trainlm', 'logsig', 'mapstd');
save(savepath, '-struct', 'networks');