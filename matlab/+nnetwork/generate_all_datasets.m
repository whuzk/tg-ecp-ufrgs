global leadnames

basedir = 'C:\physiobank\database\edb\extracted\'; 
savedir = 'C:\physiobank\database\';
ratio = 0.168;

%%
recordmap = containers.Map(leadnames, {{} {} {} ...
    {'e0207'} {} {'e0109' 'e0121' 'e0609' 'e0613'} {} {} {} {'e0403'} ...
    {'e0415' 'e0603'} {} {'e0119' 'e0161'} {'e0213' 'e0303' 'e0405'} {}});
countmap = containers.Map(leadnames, ...
    {0 0 1465 56990 0 133477 0 0 0 30548 35110 14487 108107 162747 0});
datasets = nnetwork.generate_datasets(basedir, 'Rocha', ...
    recordmap, true, countmap, [ratio 1/3]);
save([savedir 'rocha_datasets.mat'], '-struct', 'datasets');
clear datasets;

%%
recordmap = containers.Map(leadnames, {{} {} {} ...
    {} {} {'e0103' 'e0105' 'e0113' 'e0119' 'e0147'} {} {} {} {} ...
    {} {} {'e0103' 'e0105' 'e0113' 'e0119' 'e0147'} {} {}});
countmap = containers.Map(leadnames, ...
    {0 0 0 0 0 18047 0 0 0 0 0 0 18047 0 0});
datasets = nnetwork.generate_datasets(basedir, 'Mohebbi', ...
    recordmap, false, countmap, ratio);
save([savedir 'mohebbi_datasets.mat'], '-struct', 'datasets');
clear datasets;

%%
recordmap = containers.Map(leadnames, ...
    {{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}});
countmap = containers.Map(leadnames, ...
    {0 0 236 236 0 236 0 0 0 236 236 236 236 236 0});
datasets = nnetwork.generate_datasets(basedir, 'Gopalak', ...
    recordmap, true, countmap, [ratio 1/3]);
save([savedir 'gopalak_datasets.mat'], '-struct', 'datasets');
clear datasets;