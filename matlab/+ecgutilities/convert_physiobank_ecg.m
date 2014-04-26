function Result = convert_physiobank_ecg(dir, record)
% Carrega as informaçoes dos arquivos de descriçao e de dados de um ECG e
% agrega essas informaçoes numa unica variavel do MATLAB

% vai ate a pasta do registro
old_path = pwd;
cd(dir);

% extrai as informaçoes do registro
Result.Leads = wfdbdesc([record '.hea']);
M = length(Result.Leads);
for i = 1:M
    N = Result.Leads(i).LengthSamples;
    k = Result.Leads(i).RecordIndex;
    Result.Leads(i).Data = read_data(record,k,N);
end
Result.Notes = read_notes(record,M);
Result.Annotations = read_annotations(record);

% volta a pasta do usuario
cd(old_path);


function Result = read_notes(record,M)
fid = fopen([record '.hea']);
lines = textscan(fid,'%s','Delimiter','\n','Headerlines',M+1);
Result = char(lines{:});
fclose(fid);

function Result = read_data(record,si,N)
Result = zeros(N,1);
step = 100000;
i = 1;
while i <= N
    last = min(N,i+step-1);
    [~,amp] = rdsamp(record, si, last, i, 3);
    Result(i:last,:) = amp;
    i = i + step;
end

function Result = read_annotations(record)
annotators = {'atr' 'qrs' 'xyz' 'ari' '16a' 'sta' 'stb' 'stc' ...
    'pu' 'man' 'q1c' 'qt1' 'q2c' 'qt2'};
for j = 1:length(annotators)
    ann = annotators{j};
    if exist([record '.' ann],'file')
        varname = matlab.lang.makeValidName(ann);
        Result.(varname) = read_ann(record,ann);
    end
end
if exist([record '.stf'],'file')
    Result.stf = importdata([record '.stf']);
end

function Result = read_ann(record,annotator)
[sample,type,subtype,channel,num,comment] = rdann(record, annotator);
Result.Timestamp = read_time(record,sample);
Result.Sample = sample;
Result.Type = type;
Result.Subtype = subtype;
Result.Channel = channel;
Result.Num = num;
Result.Comment = cell(size(comment));
for i = 1:size(comment,1)
    Result.Comment{i} = char(comment{i});
end

function Result = read_time(record,samples)
N = length(samples);
Result = cell(N,1);
step = 3000;
i = 1;
while i <= N
    last = min(N,i+step-1);
    [ts,~] = wfdbtime(record, samples(i:last));
    Result(i:last) = ts;
    i = i + step;
end