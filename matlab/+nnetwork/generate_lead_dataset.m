function Result = generate_lead_dataset(basedir, leadname, methodname, ...
    records, discard, selcount, ratio)

files = dir([basedir '*.mat']);
list = cell(1,length(files));
rowcount = 0;
colcount = 0;
for i = 1:length(files)
    file = files(i);
    if isempty(file.name) || file.isdir
        continue;
    end
    [~,name,~] = fileparts(file.name);
    if discard == ismember(name, records)
        continue;
    end
    filepath = [basedir file.name];
    if isempty(who('-file', filepath, leadname))
        continue;
    end
    info = load(filepath, leadname);
    list{i} = info.(leadname).Datasets.(methodname);
    rowcount = rowcount + size(list{i},1);
    if colcount == 0
        colcount = size(list{i},2);
    end
end

dataset = zeros(rowcount,colcount);
iend = 0;
for i = 1:length(list)
    ibegin = iend + 1;
    iend = iend + size(list{i},1);
    dataset(ibegin:iend,:) = list{i};
end

indices = utils.select_rows(dataset, selcount, ratio);
Result = dataset(indices,:);