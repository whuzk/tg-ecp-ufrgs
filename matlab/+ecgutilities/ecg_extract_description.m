function [Result,N,M] = ecg_extract_description(lines)
%   Extrai informaçoes de um registro de ECG a partir das linhas do seu
%   arquivo de descriçao.

% cria a estrutura do ecg
Result = struct;

% faz o parse da primeira linha
info = textscan(lines{1}, '%s%s%s%s%s%s');
M = str2double(info{2}{1});
f_str = info{3}{1};
pos = strfind(f_str, '/');
if ~isempty(pos)
    f = str2double(f_str(1:pos-1));
else
    f = str2double(f_str);
end
N = str2double(info{4}{1});
len_str = datestr(datenum(0,0,0,0,0,N/f),'HH:MM:SS.FFF');
if (~isempty(info{5}) && ~isempty(info{6}))
    start_str = [info{5}{1} ' ' info{6}{1}];
elseif (~isempty(info{5}))
    start_str = info{5}{1};
else
    start_str = 'not specified';
end

% preenche as informaçoes de cabeçalho
Result.RecordName = info{1}{1};
Result.StartingTime = start_str;
Result.SignalCount = M;
Result.Signals = cell(M,1);
Result.Length = sprintf('%s (%d sample intervals)', len_str, N);
Result.SamplingFrequency = sprintf('%d Hz', f);
Result.Notes = lines(M+2:end);

% preenche as informaçoes de cada sinal
for i = 1:M
    info = textscan(lines{i+1}, '%s%s%s%s%s%s%s%s%s');
    signal.Filename = info{1}{1};
    signal.StorageFormat = info{2}{1};
    if strfind(info{3}{1}, '/')
        pos = strfind(info{3}{1}, '/');
        info{3}{1} = info{3}{1}(1:pos-1);
    end
    signal.Gain = [info{3}{1} ' adu/mV'];
    signal.Resolution = [info{4}{1} ' bits'];
    signal.ADCzero = info{5}{1};
    signal.InitialValue = info{6}{1};
    signal.Checksum = info{7}{1};
    signal.Baseline = info{8}{1};
    signal.Description = info{9}{1};
    Result.Signals{i} = signal;
    clear signal;
end
