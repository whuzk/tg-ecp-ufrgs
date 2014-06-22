function Result = get_lead_dataset(basedir, leadname, methodname, neglect)

files = dir([basedir '*.mat']);
list = cell(1,length(files));
rowcount = 0;
colcount = 0;
for i = 1:length(files)
    file = files(i);
    if isempty(file.name) || file.isdir
        continue;
    end
    filepath = [basedir file.name];
    [~,name,~] = fileparts(file.name);
    if ismember(name, neglect) || isempty(who('-file', filepath, leadname))
        continue;
    end
    info = load(filepath, leadname);
    list{i} = info.(leadname).Datasets.(methodname);
    rowcount = rowcount + size(list{i},1);
    if colcount == 0
        colcount = size(list{i},2);
    end
end

Result = zeros(rowcount,colcount);
if rowcount > 0
    iend = 0;
    for i = 1:length(list)
        ibegin = iend + 1;
        iend = iend + size(list{i},1);
        Result(ibegin:iend,:) = list{i};
    end
end