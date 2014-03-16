function Result = convert_ecg(desc_filepath, data_filepath, atr_filepath, pu_filepath)
% Carrega as informa�oes dos arquivos de descri�ao e de dados de um ECG e
% agrega essas informa�oes numa unica variavel do MATLAB
import ecgutilities.*

% verifica se existe anota�ao para o registro
atr_text = fileread(atr_filepath);
if isempty(atr_text)
    disp([atr_filepath ' : atr empty']);
    Result = [];
    return;
end

% extrai informa�oes basicas do registro
ecg_desc = importdata(desc_filepath, '\n');
[Result,~,M] = extract_description(ecg_desc);

% carrega as amostras dos sinais de ecg
ecg_data = csvread(data_filepath);
for i = 1:M
    Result.Signals{i}.Data = ecg_data(:,i+1);
end

% agrega as anota�oes do registro
Result.Annotations = extract_annotations(atr_text);
if ~isempty(pu_filepath)
    Result.PUs = cell(Result.SignalCount, 1);
    for i = 0:Result.SignalCount-1
        pu_text = fileread([pu_filepath num2str(i)]);
        Result.PUs{i+1} = extract_annotations(pu_text);
    end
end
