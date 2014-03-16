function convert_database(desc_dir, data_dir, atr_dir, ...
    pu_dir, sta_dir, stb_dir, stc_dir, save_filepath)
% carrega todos arquivos da base e converte para o formato MATLAB
import ecgutilities.*

tic;
Database = struct();

% obtem uma lista de arquivos '.hea'
desc_files = dir([desc_dir '*.hea']);

% percorre a lista e carrega os arquivos
for i = 1:length(desc_files)
    file = desc_files(i);
    [~,var_name,~] = fileparts(file.name);
    fprintf('Converting %s...\n', var_name);
    desc_filepath = [desc_dir file.name];
    data_filepath = [data_dir var_name '.sgnl'];
    atr_filepath = [atr_dir var_name '.atr'];
    if exist(pu_dir, 'dir')
        pu_filepath = [pu_dir var_name '.pu'];
    else
        pu_filepath = '';
    end
    
    ECG = convert_ecg(...
        desc_filepath, data_filepath, atr_filepath, pu_filepath);
    
    if ~isempty(ECG)
        if ~isempty(sta_dir)
            sta_filepath = [sta_dir var_name '.sta'];
            stb_filepath = [stb_dir var_name '.stb'];
            stc_filepath = [stc_dir var_name '.stc'];
            ECG.ST = convert_st_annot(...
                sta_filepath, stb_filepath, stc_filepath);
        end
        Database.(var_name) = ECG;
    end
end

% salva toda a base para o arquivo de saida
fprintf('Saving data to file...\n');
save(save_filepath, '-struct', 'Database');
fprintf('Saved to ''%s''\n', save_filepath);

toc;