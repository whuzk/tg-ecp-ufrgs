function [Fs,Sc,Leads] = interpret(Record)
% Fs - sampling frequency
% Sc - number of leads
% Leads - the signals and their main properties

Fs = sscanf(Record.SamplingFrequency,'%d');
idx = isqrs(char(Record.Annotations.Type));
qrs = str2double(Record.Annotations.Sample(idx));
Sc = Record.SignalCount;
Leads = cell(1,Sc);
for i = 1:Sc
    Lead = Record.Signals{i};
    Leads{i}.name = Lead.Description;
    Leads{i}.bres = sscanf(Lead.Resolution,'%d');
    Leads{i}.gain = sscanf(Lead.Gain,'%d');
    Leads{i}.zero = sscanf(Lead.ADCzero,'%d');
    Leads{i}.qrsi = qrs;
    if max(Lead.Data) > 5
        Leads{i}.data = Lead.Data;
        Leads{i}.iniv = Lead.Data(1);
    else
        Leads{i}.data = round(Lead.Data * Leads{i}.gain);
        Leads{i}.iniv = round(Lead.Data(1) * Leads{i}.gain);
    end
end

function Result = isqrs(type)
Result = ismember(type, 'NLRBAaJSVFejnE/fQ?');