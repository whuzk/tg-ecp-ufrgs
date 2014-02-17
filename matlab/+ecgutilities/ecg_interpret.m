function [Fs,Bp,Sc,Leads,Data] = ecg_interpret(ECG)

Fs = sscanf(ECG.SamplingFrequency,'%d');
Bp = str2double(ECG.Annotations.Sample);
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
