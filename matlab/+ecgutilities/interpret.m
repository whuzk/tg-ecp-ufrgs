function [Fs,Bp,Sc,Leads] = interpret(ECG)
% Fs - sampling frequency
% Bp - beat locations
% Sc - number of leads
% Leads - the signals

% get sampling frequency
Fs = sscanf(ECG.SamplingFrequency,'%d');

% get beat locations
qrs = isqrs(char(ECG.Annotations.Type));
Bp = str2double(ECG.Annotations.Sample(qrs));

% get the signals
Sc = ECG.SignalCount;
Leads = cell(1,Sc);
for i = 1:Sc
    Lead = ECG.Signals{i};
    Leads{i}.name = Lead.Description;
    Leads{i}.res = sscanf(Lead.Resolution,'%d');
    Leads{i}.gain = sscanf(Lead.Gain,'%d');
    Leads{i}.data = Lead.Data;
end

function Result = isqrs(type)
Result = ismember(type, 'NLRBAaJSVFejnE/fQ?');