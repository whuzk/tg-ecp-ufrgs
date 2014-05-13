function group_physiobank_db(db_dir, save_path)
import ecgutilities.*

% obtem uma lista de arquivos '.mat'
desc_files = dir([db_dir filesep '*.mat']);

% percorre a lista e carrega os arquivos
for i = 1:length(desc_files)
    file = desc_files(i);
    [~,record,~] = fileparts(file.name);
    fprintf('Loading %s...\n', record);
    var_name = matlab.lang.makeValidName(record);
    database.(var_name) = load([db_dir filesep file.name]);
end

% salva em arquivo
fprintf('Saving data to file...\n');
save(save_path, '-struct', 'database');
fprintf('Saved to ''%s''\n', save_path);