function [Fs,Bp,Sc,Leads,Data] = interpret(ECG)

Fs = sscanf(ECG.SamplingFrequency,'%d');
index = ...
    ~strcmp(ECG.Annotations.Type, '+') &...
    ~strcmp(ECG.Annotations.Type, '|') &...
    ~strcmp(ECG.Annotations.Type, '~') &...
    ~strcmp(ECG.Annotations.Type, '"') &...
    ~strcmp(ECG.Annotations.Type, 's') &...
    ~strcmp(ECG.Annotations.Type, 'x') &...
    ~strcmp(ECG.Annotations.Type, 'T');
Bp = str2double(ECG.Annotations.Sample(index));
N = sscanf(ECG.Length(strfind(ECG.Length,'(')+1:end),'%d');
Sc = ECG.SignalCount;

Leads = cell(1,Sc);
Data = zeros(N,Sc);
for i = 1:Sc
    Lead = ECG.Signals{i};
    Leads{i} = Lead.Description;
    Gain = sscanf(Lead.Gain,'%d');
    Signal = Lead.Data-Lead.Data(1);
    Span = (max(Signal)-min(Signal))/2;
    if (Span <= 24)
        Gain = 1;
    end
    Data(:,i) = Signal/Gain;
end
