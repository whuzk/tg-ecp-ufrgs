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

indices = select_rows(dataset, selcount, ratio);
Result = dataset(indices,:);


function Result = select_rows(dataset, selcount, ratio)
if isempty(ratio)
    if size(dataset,1) < selcount
        error('number of detected beats is less than expected: %d', rowcount);
    end
    Result = sort(randperm(size(dataset,1), selcount));
    
elseif length(ratio) == 1
    normalcount = ceil((1-ratio)*selcount);
    ischemiccount = floor(ratio*selcount);
    
    normal = find(dataset(:,end-1) == 0 & dataset(:,end) == 0);
    ischemic = find(dataset(:,end-1) ~= 0 | dataset(:,end) ~= 0);
    
    if length(normal) < normalcount
        error('number of normal beats is less than expected: %d', length(normal));
    elseif length(ischemic) < ischemiccount
        error('number of ischemic beats is less than expected: %d', length(ischemic));
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(ischemic), ischemiccount);
    Result = unique([normal(indices1); ischemic(indices2)]);
    
elseif length(ratio) == 2
    normalcount = ceil((1-ratio(1))*selcount);
    ischemiccount = floor(ratio(1)*selcount);
    sttcount = floor(ratio(2)*ischemiccount);
    stcount = ceil((ischemiccount-sttcount)/2);
    tcount = ischemiccount - sttcount - stcount;
    
    normal = find(dataset(:,end-1) == 0 & dataset(:,end) == 0);
    stchange = find(dataset(:,end-1) ~= 0 & dataset(:,end) == 0);
    tchange = find(dataset(:,end-1) == 0 & dataset(:,end) ~= 0);
    sttchange = find(dataset(:,end-1) ~= 0 & dataset(:,end) ~= 0);
    
    if length(normal) < normalcount
        error('number of normal beats is less than expected: %d', length(normal));
    elseif length(stchange) < stcount
        error('number of beats with ST change is less than expected: %d', length(stchange));
    elseif length(tchange) < tcount
        error('number of beats with T change is less than expected: %d', length(tchange));
    elseif length(sttchange) < sttcount
        error('number of beats with ST-T change is less than expected: %d', length(sttchange));
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(stchange), stcount);
    indices3 = randperm(length(tchange), tcount);
    indices4 = randperm(length(sttchange), sttcount);
    
    Result = unique([normal(indices1); stchange(indices2); ...
        tchange(indices3); sttchange(indices4)]);
    
else
    error('ratio must be a vector of length <= 2');
end