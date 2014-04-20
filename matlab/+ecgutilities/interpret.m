function [Fs,Bp,Sc,Leads] = interpret(ECG)
% Fs - sampling frequency
% Bp - beat locations
% Sc - number of leads
% Leads - the signals

% get basic info
Fs = sscanf(ECG.SamplingFrequency,'%d');
N = sscanf(ECG.Length(strfind(ECG.Length,'(')+1:end),'%d');
Sc = ECG.SignalCount;

% non-beat annotation characters
exc = {
    '['     % Start of ventricular flutter/fibrillation
    '!'     % Ventricular flutter wave
    ']'     % End of ventricular flutter/fibrillation
    'x'     % Non-conducted P-wave (blocked APC)
    '('     % Waveform onset
    ')'     % Waveform end
    'p'     % Peak of P-wave
    't'     % Peak of T-wave
    'u'     % Peak of U-wave
    '`'     % PQ junction
    ''''    % J-point
    '^'     % (Non-captured) pacemaker artifact
    '|'     % Isolated QRS-like artifact
    '~'     % Change in signal quality
    '+'     % Rhythm change
    's'     % ST segment change
    'T'     % T-wave change
    '*'     % Systole
    'D'     % Diastole
    '='     % Measurement annotation
    '"'     % Comment annotation
    '@'     % Link to external data
};

% get beat locations
index = true(length(ECG.Annotations.Type),1);
for i = 1:length(exc)
    index = index & ~strcmp(ECG.Annotations.Type, exc{i});
end
Bp = str2double(ECG.Annotations.Sample(index));

% get the signals
Leads = cell(1,Sc);
for i = 1:Sc
    Lead = ECG.Signals{i};
    Leads{i}.name = Lead.Description;
    Leads{i}.res = sscanf(Lead.Resolution,'%d');
    Leads{i}.gain = sscanf(Lead.Gain,'%d');
    Leads{i}.data = Lead.Data;
end
