function convert_physiobank_db(source_dir, target_dir)
% carrega todos arquivos da base e converte para o formato MATLAB
import ecgutilities.*

tic;

% cria o diretorio destino, se necessario
if ~exist(target_dir,'dir') && ~mkdir(target_dir)
    error('could not create new directory: %s', target_dir);
end

% obtem uma lista de arquivos '.hea'
desc_files = dir([source_dir filesep '*.hea']);

% percorre a lista e carrega os arquivos
for i = 1:length(desc_files)
    file = desc_files(i);
    [~,record,~] = fileparts(file.name);
    
    % converte o ecg
    fprintf('Converting %s...\n', record);
    var_name = matlab.lang.makeValidName(record);
    dum.(var_name) = convert_physiobank_ecg(source_dir, record);
    
    % salva em arquivo
    save_file = [target_dir filesep record '.mat'];
    fprintf('Saving data to file...\n');
    save(save_file,'-struct','dum');
    fprintf('Saved to ''%s''\n', save_file);
end

toc;