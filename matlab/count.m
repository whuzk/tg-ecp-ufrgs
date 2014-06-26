basedir = 'C:\physiobank\database\edb\extracted\'; 
files = dir([basedir '*.mat']);
list = cell(1,length(files));
ischcount = 0;
totalcount = 0;
for i = 1:length(files)
    file = files(i);
    if isempty(file.name) || file.isdir
        continue;
    end
    info = load([basedir file.name]);
    fields = fieldnames(info);
    dataset = info.(fields{2}).Datasets.Rocha;
    ischcount = ischcount + length(find(dataset(:,end-1)~=0 | dataset(:,end)~=0));
    totalcount = totalcount + size(dataset,1);
end
ischcount/totalcount