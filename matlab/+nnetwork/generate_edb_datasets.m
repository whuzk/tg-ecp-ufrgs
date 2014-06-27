global edbleadnames

basedir = 'C:\physiobank\database\edb\extracted\'; 
savedir = 'C:\physiobank\database\';

stratiomap = containers.Map(edbleadnames, {...
    [0.9262    0.0738         0]
    [0.8943    0.0339    0.0718]
    [0.9307    0.0337    0.0356]
    [0.8463    0.0044    0.1493]
    [0.8608    0.0175    0.1217]
    [0.9492    0.0012    0.0496]
    [0.8723    0.0082    0.1195]
    [0.8484    0.0061    0.1455]
});

tratiomap = containers.Map(edbleadnames, {...
    [0.8341    0.1659         0]
    [0.9690    0.0074    0.0236]
    [0.9264    0.0429    0.0307]
    [0.8868    0.0338    0.0795]
    [0.8050    0.0236    0.1714]
    [0.8925    0.0713    0.0362]
    [0.8868    0.0532    0.0600]
    [0.9491    0.0240    0.0269]
});

recordmap = containers.Map(edbleadnames, {...
    {}
    {}
    {}
    {}
    {}
    {}
    {}
    {}
});

countmap = containers.Map(edbleadnames, {...
    5000
    50000
    100000
    30000
    30000
    10000
    100000
    100000
});

Rocha = nnetwork.generate_datasets(edbleadnames, ...
    basedir, 'Rocha', recordmap, true, countmap, stratiomap, tratiomap);
save([savedir 'rocha_datasets.mat'], '-struct', 'Rocha');

Mohebbi = nnetwork.generate_datasets(edbleadnames, ...
    basedir, 'Mohebbi', recordmap, true, countmap, stratiomap, tratiomap);
save([savedir 'mohebbi_datasets.mat'], '-struct', 'Mohebbi');

Gopalak = nnetwork.generate_datasets(edbleadnames, ...
    basedir, 'Gopalak', recordmap, true, countmap, stratiomap, tratiomap);
save([savedir 'gopalak_datasets.mat'], '-struct', 'Gopalak');