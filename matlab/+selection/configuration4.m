import selection.*
global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\';
savepath = 'C:\physiobank\database\superset4.mat';
load('matfiles\configuration4.mat');

RochaSet = generate_datasets(edbleadnames, basedir, 'Rocha', ...
    strecordmap, false, countmap, stratiomap, 0);
MohebbiSet = generate_datasets(edbleadnames, basedir, 'Mohebbi', ...
    strecordmap, false, countmap, stratiomap, 0);
GopalakSTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    strecordmap, false, countmap, stratiomap, 0);
GopalakTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    trecordmap, false, countmap, tratiomap, 1);

save(savepath, 'RochaSet', 'MohebbiSet', 'GopalakSTSet', 'GopalakTSet');
clear RochaSet MohebbiSet GopalakSTSet GopalakTSet