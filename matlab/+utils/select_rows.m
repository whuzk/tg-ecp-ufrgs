function Result = select_rows(dataset, selcount, ratio, twave)

if length(ratio) == 1
    if size(dataset,1) < selcount
        error('number of detected beats is less than expected: %d < %d',...
            size(dataset,1), selcount);
    end
    Result = sort(randperm(size(dataset,1), selcount));
    
elseif length(ratio) == 2
    normalcount = ceil(ratio(1)*selcount);
    ischemiccount = selcount - normalcount;
    
    diagidx = size(dataset,2) - double(~twave);
    normal = find(dataset(:,diagidx) == 0);
    ischemic = find(dataset(:,diagidx) ~= 0);
    
    if length(normal) < normalcount
        error('number of normal beats is less than expected: %d < %d', ...
            length(normal), normalcount);
    elseif length(ischemic) < ischemiccount
        error('number of ischemic beats is less than expected: %d < %d', ...
            length(ischemic), ischemiccount);
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(ischemic), ischemiccount);
    Result = unique([normal(indices1); ischemic(indices2)]);
    
elseif length(ratio) == 3
    normalcount = ceil(ratio(1)*selcount);
    ischemiccount = selcount - normalcount;
    elevcount = ceil(ratio(2)/(ratio(2)+ratio(3))*ischemiccount);
    depcount = ischemiccount - elevcount;
    
    diagidx = size(dataset,2) - double(~twave);
    normal = find(dataset(:,diagidx) == 0);
    elev = find(dataset(:,diagidx) > 0);
    dep = find(dataset(:,diagidx) < 0);
    
    if length(normal) < normalcount
        error('number of normal beats is less than expected: %d < %d', ...
            length(normal), normalcount);
    elseif length(elev) < elevcount
        error('number of beats with elevation is less than expected: %d < %d', ...
            length(elev), elevcount);
    elseif length(dep) < depcount
        error('number of beats with depression is less than expected: %d < %d', ...
            length(dep), depcount);
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(elev), elevcount);
    indices3 = randperm(length(dep), depcount);
    Result = unique([normal(indices1); elev(indices2); dep(indices3)]);
    
else
    error('length of ratio must be between 1 and 3');
end