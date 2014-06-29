import nnetwork.*
global edbleadnames

savepath = 'C:\physiobank\database\networks3.mat';
load('C:\physiobank\database\superset3.mat');

RochaNet = train_networks(edbleadnames, ...
    RochaSet, 'trainlm', 'tansig', 'mapstd');
MohebbiNet = train_networks(edbleadnames, ...
    MohebbiSet, 'trainlm', 'tansig', 'mapstd');
GopalakSTNet = train_networks(edbleadnames, ...
    GopalakSTSet, 'trainlm', 'tansig', 'mapstd');
GopalakTNet = train_networks(edbleadnames, ...
    GopalakTSet, 'trainlm', 'tansig', 'mapstd');

save(savepath, 'RochaNet', 'MohebbiNet', 'GopalakSTNet', 'GopalakTNet');