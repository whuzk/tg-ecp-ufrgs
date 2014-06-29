import selection.*
global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\';
savepath = 'C:\physiobank\database\superset1.mat';
load('matfiles\configuration1.mat');

RochaSet = generate_datasets(edbleadnames, basedir, 'Rocha', ...
    rocharecordmap, true, rochacountmap, stratiomap, 0);
MohebbiSet = generate_datasets(edbleadnames, basedir, 'Mohebbi', ...
    mohebbirecordmap, false, mohebbicountmap, mohebbiratiomap, 0);
GopalakSTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    gopalakrecordmap, true, gopalakcountmap, stratiomap, 0);
GopalakTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    gopalakrecordmap, true, gopalakcountmap, tratiomap, 1);

save(savepath, 'RochaSet', 'MohebbiSet', 'GopalakSTSet', 'GopalakTSet');