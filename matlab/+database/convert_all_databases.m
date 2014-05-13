function convert_all_databases
% converte as bases de dados de ECGs

source_dir = 'I:\AppData\physiobank\database';
target_dir = 'C:\physiobank\database';
selected = {'edb', 'mitdb', 'qtdb'};

dirs = dir(source_dir);
for i = 1:length(dirs)
    file = dirs(i);
    if file.isdir && ismember(file.name, selected)
        tic;
        fprintf('Converting %s database...\n', file.name);
        source = [source_dir filesep file.name];
        target = [target_dir filesep file.name];
        database.convert_physiobank_db(source, target);
        fprintf('Saved to %s\n', target);
        toc;
    end
end