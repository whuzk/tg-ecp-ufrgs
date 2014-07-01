import selection.*
global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\';
savepath = 'C:\physiobank\database\superset3.mat';
load('matfiles\configuration3.mat');

RochaSet = generate_datasets(edbleadnames, basedir, 'Rocha', ...
    strecordmap, false, countmap, ratiomap, 0);
MohebbiSet = generate_datasets(edbleadnames, basedir, 'Mohebbi', ...
    strecordmap, false, countmap, ratiomap, 0);
GopalakSTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    strecordmap, false, countmap, ratiomap, 0);
GopalakTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    trecordmap, false, countmap, ratiomap, 1);

save(savepath, 'RochaSet', 'MohebbiSet', 'GopalakSTSet', 'GopalakTSet');
clear RochaSet MohebbiSet GopalakSTSet GopalakTSet