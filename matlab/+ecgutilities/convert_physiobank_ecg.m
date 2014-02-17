function Result = convert_physiobank_ecg(db_dir, record_name)
% Carrega as informaçoes dos arquivos de descriçao e de dados de um ECG e
% agrega essas informaçoes numa unica variavel do MATLAB

% extrai informaçoes basicas do registro
ecg_desc = importdata([db_dir record_name '.hea'], '\n');
[Result,N,M] = utilities.ecg_extract_description(ecg_desc);

% vai ate a pasta do registro
old_path = pwd;
cd(db_dir);

% carrega as amostras dos sinais de ecg
ecg_data = zeros(N,M);
step = 100000;
a = floor(N/step);
rest = N - a*step;
for i = 1:a
    i0 = step*(i-1)+1;
    i1 = step*i;
    [~,ecg_data(i0:i1,1:M)] = rdsamp(record_name, 1:M, i1, i0);
end
if rest > 0
    i0 = a*step+1;
    [~,ecg_data(i0:N,1:M)] = rdsamp(record_name, 1:M, N, i0);
end
for i = 1:M
    Result.Signals{i}.Data = ecg_data(:,i);
end

% agrega as anotaçoes do ECG
[ann,type,subtype,chan,num,comments] = rdann(record_name, 'atr');
Result.Annotations.DateTime = wfdbtime(record_name,ann);
Result.Annotations.Sample = strtrim(cellstr(num2str(ann)));
Result.Annotations.Type = cellstr(type);
Result.Annotations.Sub = cellstr(subtype);
Result.Annotations.Chan = strtrim(cellstr(num2str(chan)));
Result.Annotations.Num = strtrim(cellstr(num2str(num)));
Result.Annotations.Aux = cell(size(comments,1),1);
for i = 1:size(comments,1)
    Result.Annotations.Aux{i} = char(comments{i});
end

% volta a pasta do usuario
cd(old_path);
