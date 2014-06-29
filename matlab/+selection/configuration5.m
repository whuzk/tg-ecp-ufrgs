import selection.*
global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\';
savepath = 'C:\physiobank\database\superset5.mat';
load('matfiles\configuration5.mat');

RochaSet = generate_datasets(edbleadnames, basedir, 'Rocha', ...
    recordmap, true, countmap, ratiomap, 0);
MohebbiSet = generate_datasets(edbleadnames, basedir, 'Mohebbi', ...
    recordmap, true, countmap, ratiomap, 0);
GopalakSTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    recordmap, true, countmap, ratiomap, 0);
GopalakTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    recordmap, true, countmap, ratiomap, 1);

save(savepath, 'RochaSet', 'MohebbiSet', 'GopalakSTSet', 'GopalakTSet');