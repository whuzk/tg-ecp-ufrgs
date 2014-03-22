function [Fs,Bp,Sc,Leads,Data] = interpret(ECG)

Fs = sscanf(ECG.SamplingFrequency,'%d');
index = ...
    ~strcmp(ECG.Annotations.Type, '+') &...
    ~strcmp(ECG.Annotations.Type, '|') &...
    ~strcmp(ECG.Annotations.Type, '~') &...
    ~strcmp(ECG.Annotations.Type, 's') &...
    ~strcmp(ECG.Annotations.Type, 'T');
Bp = str2double(ECG.Annotations.Sample(index));
N = sscanf(ECG.Length(strfind(ECG.Length,'(')+1:end),'%d');
Sc = ECG.SignalCount;

Leads = cell(1,Sc);
Data = zeros(N,Sc);
for i = 1:Sc
    Signal = ECG.Signals{i};
    Leads{i} = Signal.Description;
    Gain = sscanf(Signal.Gain,'%d');
    Offset = sscanf(Signal.InitialValue,'%d');
    Data(:,i) = (Signal.Data - Offset) / Gain;
end
