import selection.*
global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\';
savepath = 'C:\physiobank\database\superset2.mat';
load('matfiles\configuration2.mat');

RochaSet = generate_datasets(edbleadnames, basedir, 'Rocha', ...
    rocharecordmap, false, countmap, ratiomap, 0);
MohebbiSet = generate_datasets(edbleadnames, basedir, 'Mohebbi', ...
    mohebbirecordmap, false, countmap, ratiomap, 0);
GopalakSTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    gopalakstrecordmap, false, countmap, ratiomap, 0);
GopalakTSet = generate_datasets(edbleadnames, basedir, 'Gopalak', ...
    gopalaktrecordmap, false, countmap, ratiomap, 1);

save(savepath, 'RochaSet', 'MohebbiSet', 'GopalakSTSet', 'GopalakTSet');
clear RochaSet MohebbiSet GopalakSTSet GopalakTSet