function Result = convert_physiobank_ecg(dir, record)
% Carrega as informaçoes dos arquivos de descriçao e de dados de um ECG e
% agrega essas informaçoes numa unica variavel do MATLAB

% vai ate a pasta do registro
old_path = pwd;
cd(dir);

% extrai as informaçoes do registro
Result.Info = wfdbdesc([record '.hea']);
Result.Notes = read_notes(record);
Result.Annotations = read_annotations(record);
N = Result.Info(1).LengthSamples;
Result.SignalCount = length(Result.Info);
Result.SignalData = read_data(record,N,Result.SignalCount);

% volta a pasta do usuario
cd(old_path);


function Result = read_notes(record)
fid = fopen([record '.hea']);
line = fgets(fid);
Result = '';
while ischar(line)
    if strfind(line,'#') == 1
        Result = [Result line(2:end)];
    end
    line = fgets(fid);
end
fclose(fid);

function Result = read_data(record,N,Sc)
Result = zeros(N,Sc);
step = 100000;
i = 1;
while i <= N
    last = min(N,i+step-1);
    [~,amp] = rdsamp(record, [], last, i, 3);
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
    Result.stf = read_stf([record '.stf']);
end

function Result = read_ann(record,annotator)
[Sample,Type,Subtype,Channel,Num,Aux] = rdann(record, annotator);
Timestamp = read_time(record,Sample);
Comment = cell(size(Aux));
for i = 1:size(Aux,1)
    Comment{i} = char(Aux{i});
end
Result = table(Timestamp,Sample,Type,Subtype,Channel,Num,Comment);

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

function Result = read_stf(filepath)
data = importdata(filepath);
Sample = data(:,1);
M = (size(data,2)-1)/3;
Result = cell(1,M);
for i = 0:M-1
    STLevelFunction = data(:,3*i+2);
    STLevelReference = data(:,3*i+3);
    STLevelDeviation = data(:,3*i+4);
    Result{i+1} = table(Sample, ...
        STLevelFunction, STLevelReference, STLevelDeviation);
end