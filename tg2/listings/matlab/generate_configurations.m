function generate_configurations
global edbleadnames

basedir = 'C:\physiobank\database\edb\';
count = utils.count_classes(basedir);
load('matfiles\indstats.mat');
lists = get_selection_lists(tables);
merge = get_merged_selection(lists);

%% configuration 1
[stratiomap,tratiomap] = get_edb_ratio_maps(count);
temp = cell2mat(stratiomap.values');
temp = [temp(:,1) temp(:,2)+temp(:,3)];
mohebbiratiomap = containers.Map(edbleadnames, num2cell(temp,2));
[rocharecordmap,mohebbirecordmap,gopalakrecordmap] = ...
    get_original_record_maps();
[rochacountmap,mohebbicountmap,gopalakcountmap] = get_original_count_maps();

save('matfiles\configuration1.mat', 'stratiomap', 'tratiomap', ...
    'mohebbiratiomap', 'rocharecordmap', 'mohebbirecordmap', ...
    'gopalakrecordmap', 'rochacountmap', 'mohebbicountmap', 'gopalakcountmap');

%% configuration 2
ratiomap = containers.Map(edbleadnames, {[] [] [] [] [] [] [] []});
[rocharecordmap,mohebbirecordmap,gopalakstrecordmap,gopalaktrecordmap] = ...
    get_proposed_record_maps(lists, basedir);
countmap = containers.Map(edbleadnames, {[] [] [] [] [] [] [] []});

save('matfiles\configuration2.mat', 'ratiomap', 'rocharecordmap', ...
    'mohebbirecordmap', 'gopalakstrecordmap', 'gopalaktrecordmap', ...
    'countmap');

%% configuration 3
[strecordmap,trecordmap] = get_merged_record_maps(merge, basedir);
countmap = containers.Map(edbleadnames, {[] [] [] [] [] [] [] []});

save('matfiles\configuration3.mat', 'ratiomap', 'strecordmap', ...
    'trecordmap', 'countmap');

%% configuration 4
[stratiomap,tratiomap] = get_edb_ratio_maps(count);
countmap = containers.Map(edbleadnames, ...
    {5000 20000 100000 30000 20000 10000 100000 100000});

save('matfiles\configuration4.mat', 'stratiomap', 'tratiomap', ...
    'strecordmap', 'trecordmap', 'countmap');

%% configuration 5
ratiomap = containers.Map(edbleadnames, {[] [] [] [] [] [] [] []});
recordmap = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
countmap = containers.Map(edbleadnames, {[] [] [] [] [] [] [] []});

save('matfiles\configuration5.mat', 'ratiomap', 'recordmap', 'countmap');

function [stratiomap,tratiomap] = get_edb_ratio_maps(count)
global edbleadnames
N = length(edbleadnames);
ST = zeros(N,3);
T = zeros(N,3);
for i = 1:N
    lead = edbleadnames{i};
    total = sum(sum(count.(lead)));
    ST(i,:) = sum(count.(lead),2)/total;
    T(i,:) = sum(count.(lead),1)/total;
end
stratiomap = containers.Map(edbleadnames, num2cell(ST,2));
tratiomap = containers.Map(edbleadnames, num2cell(T,2));

function [rocha,mohebbi,gopalak] = get_original_record_maps()
global edbleadnames
rocha = containers.Map(edbleadnames, {...
    {}
    {'e0207'}
    {'e0109' 'e0121' 'e0609' 'e0613'}
    {'e0403'}
    {'e0415' 'e0603'}
    {}
    {'e0119' 'e0161'}
    {'e0207' 'e0213' 'e0303' 'e0405'}
});
mohebbi = containers.Map(edbleadnames, {...
    {}
    {}
    {'e0103' 'e0105' 'e0113' 'e0119' 'e0147'}
    {}
    {}
    {}
    {'e0103' 'e0105' 'e0113' 'e0119' 'e0147'}
    {}
});
gopalak = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});

function [rocha,mohebbi,gopalak] = get_original_count_maps()
global edbleadnames
rocha = containers.Map(edbleadnames, ...
    {1465 56990 133477 30548 35110 14487 108107 162747});
mohebbi = containers.Map(edbleadnames, ...
    {0 0 18047 0 0 0 18047 0});
gopalak = containers.Map(edbleadnames, ...
    {236 236 236 236 236 236 236 236});

function lists = get_selection_lists(tables)
lists = cell(length(tables),1);
for i = 1:length(tables)
    mytable = tables{i};
    sel1 = mytable.Sensitivity >= 0.8;
    sel2 = mytable.PositivePred >= 0.8;
    lists{i} = mytable(sel1 & sel2, 'RecordName');
end

function merge = get_merged_selection(lists)
merge = cell(4,1);
for i = 1:2
    list1 = lists{3*(i-1)+1};
    list2 = lists{3*(i-1)+2};
    list3 = lists{3*(i-1)+3};
    list4 = union(union(list1{:,:},list2{:,:}),list3{:,:});
    idx = true(size(list4));
    for j = 1:length(list4)
        name = list4{j};
        cond1 = any(strcmp(name, list1{:,:}));
        cond2 = any(strcmp(name, list2{:,:}));
        cond3 = any(strcmp(name, list3{:,:}));
        if cond1 + cond2 + cond3 <= 1
            idx(j) = false;
        end
    end
    merge{i} = array2table(list4(idx),'VariableNames',{'RecordName'});
end
merge{3} = lists{end-1};
merge{4} = lists{end};

function [rocha,mohebbi,gopalakst,gopalakt] = get_proposed_record_maps(lists, basedir)
global edbleadnames
rocha = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
mohebbi = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
gopalakst = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
gopalakt = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
for i = 1:2
    list = lists{3*(i-1)+1}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        rocha(lead) = [rocha(lead); list{j}];
    end
    list = lists{3*(i-1)+2}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        mohebbi(lead) = [mohebbi(lead); list{j}];
    end
    list = lists{3*(i-1)+3}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        gopalakst(lead) = [gopalakst(lead); list{j}];
    end
end
for i = 1:2
    list = lists{i+6}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        gopalakt(lead) = [gopalakt(lead); list{j}];
    end
end

function [strecordmap,trecordmap] = get_merged_record_maps(merge, basedir)
global edbleadnames
strecordmap = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
trecordmap = containers.Map(edbleadnames, {{} {} {} {} {} {} {} {}});
for i = 1:2
    list = merge{i}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        strecordmap(lead) = [strecordmap(lead); list{j}];
    end
    list = merge{i+2}{:,:};
    for j = 1:length(list)
        record = load([basedir list{j} '.mat']);
        lead = record.Info(i).Description;
        trecordmap(lead) = [trecordmap(lead); list{j}];
    end
end