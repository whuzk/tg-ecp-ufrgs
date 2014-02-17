function Result = convert_st_annot(sta_filepath, stb_filepath, stc_filepath)
% Converte as anotaçoes de desvio ST do registro

% carrega os arquivos de anotaçoes
sta_text = fileread(sta_filepath);
stb_text = fileread(stb_filepath);
stc_text = fileread(stc_filepath);

% verifica se existe anotaçao para o registro
if isempty(sta_text) || isempty(stb_text) || isempty(stc_text)
    disp([atr_filepath ' : atr empty']);
    Result = [];
    return;
end

% converte as anotaçoes
Result.sta = utilities.ecg_get_annotations(sta_text);
Result.stb = utilities.ecg_get_annotations(stb_text);
Result.stc = utilities.ecg_get_annotations(stc_text);
