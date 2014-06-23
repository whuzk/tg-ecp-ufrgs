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