% main.m
%   Rotina principal do trabalho.
%

% carrega e processa todos os registros da base EDB
basedir = 'C:\physiobank\database\edb\';
mkdir([basedir 'extracted']);
files = dir([basedir '*.mat']);
for i = 1:length(files)
    file = files(i);
    if isempty(file.name) || file.isdir
        continue;
    end
    [~,name,~] = fileparts(file.name);
    fprintf('Processing %s...\n', name);
    record = load([basedir file.name]);
    var = struct;
    for j = 1:record.SignalCount
        [info,leadname] = utils.extract_datasets(record, j);
        leadname = matlab.lang.makeValidName(leadname);
        var.(leadname) = info;
    end
    filename = [basedir 'extracted\' file.name];
    fprintf('Saving %s...\n', name);
    save(filename, '-struct', 'var');
    clear var;
end