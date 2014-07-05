function Result = select_rows(dataset, selcount, ratio, feature)

if length(ratio) == 1
    if size(dataset,1) < selcount
        warning('number of beats is less than expected: %d < %d', ...
            size(dataset,1), selcount);
        selcount = size(dataset,1);
    end
    Result = sort(randperm(size(dataset,1), selcount));
    
elseif length(ratio) == 2
    normalcount = ceil(ratio(1)*selcount);
    ischemiccount = selcount - normalcount;
    
    normal = find(dataset(:,end-1+feature) == 0);
    ischemic = find(dataset(:,end-1+feature) ~= 0);
    
    if length(normal) < normalcount
        warning('number of normal beats is less than expected: %d < %d', ...
            length(normal), normalcount);
        normalcount = length(normal);
        ischemiccount = ceil(normalcount*ratio(2)/ratio(1));
        fprintf('new count: %d\n', normalcount + ischemiccount);
    end
    if length(ischemic) < ischemiccount
        warning('number of ischemic beats is less than expected: %d < %d', ...
            length(ischemic), ischemiccount);
        ischemiccount = length(ischemic);
        normalcount = ceil(ischemiccount*ratio(1)/ratio(2));
        fprintf('new count: %d\n', normalcount + ischemiccount);
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(ischemic), ischemiccount);
    Result = unique([normal(indices1); ischemic(indices2)]);
    
elseif length(ratio) == 3
    normalcount = ceil(ratio(1)*selcount);
    ischemiccount = selcount - normalcount;
    elevcount = ceil(ratio(2)/(ratio(2)+ratio(3))*ischemiccount);
    depcount = ischemiccount - elevcount;
    
    normal = find(dataset(:,end-1+feature) == 0);
    elev = find(dataset(:,end-1+feature) > 0);
    dep = find(dataset(:,end-1+feature) < 0);
    
    if length(normal) < normalcount
        warning('number of normal beats is less than expected: %d < %d', ...
            length(normal), normalcount);
        normalcount = length(normal);
        elevcount = ceil(normalcount*ratio(2)/ratio(1));
        depcount = ceil(normalcount*ratio(3)/ratio(1));
        fprintf('new count: %d\n', normalcount + elevcount + depcount);
    end
    if length(elev) < elevcount
        warning('number of beats with elevation is less than expected: %d < %d', ...
            length(elev), elevcount);
        elevcount = length(elev);
        normalcount = ceil(elevcount*ratio(1)/ratio(2));
        depcount = ceil(elevcount*ratio(3)/ratio(2));
        fprintf('new count: %d\n', normalcount + elevcount + depcount);
    end
    if length(dep) < depcount
        warning('number of beats with depression is less than expected: %d < %d', ...
            length(dep), depcount);
        depcount = length(dep);
        normalcount = ceil(depcount*ratio(1)/ratio(3));
        elevcount = ceil(depcount*ratio(2)/ratio(3));
        fprintf('new count: %d\n', normalcount + elevcount + depcount);
    end
    
    indices1 = randperm(length(normal), normalcount);
    indices2 = randperm(length(elev), elevcount);
    indices3 = randperm(length(dep), depcount);
    Result = unique([normal(indices1); elev(indices2); dep(indices3)]);
    
else
    error('length of ratio must be between 1 and 3');
end