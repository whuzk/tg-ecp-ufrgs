import nnetwork.*;
basedir = 'C:\physiobank\database\';

%%
loadpath = [basedir 'rocha_datasets.mat'];
savepath = [basedir 'rocha_networks.mat'];
train_networks(loadpath, savepath, 'trainlm');

%%
loadpath = [basedir 'mohebbi_datasets.mat'];
savepath = [basedir 'mohebbi_networks.mat'];
train_networks(loadpath, savepath, 'trainlm');

%%
loadpath = [basedir 'gopalak_datasets.mat'];
savepath = [basedir 'gopalak_networks.mat'];
train_networks(loadpath, savepath, 'trainlm');