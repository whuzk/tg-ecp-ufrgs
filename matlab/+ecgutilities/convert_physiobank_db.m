function convert_physiobank_db(db_dir, save_filepath)
% carrega todos arquivos da base e converte para o formato MATLAB
import ecgutilities.*

tic;
Database = struct();

% obtem uma lista de arquivos '.hea'
desc_files = dir([db_dir '*.hea']);

% percorre a lista e carrega os arquivos
for i = 1:length(desc_files)
    file = desc_files(i);
    [~,var_name,~] = fileparts(file.name);
    fprintf('Converting %s...\n', var_name);
    ECG = convert_physiobank_ecg(db_dir, var_name);
    var_name = matlab.lang.makeValidName(var_name);
    Database.(var_name) = ECG;
end

% salva toda a base para o arquivo de saida
fprintf('Saving data to file...\n');
save(save_filepath, '-struct', 'Database');
fprintf('Saved to ''%s''\n', save_filepath);

toc;