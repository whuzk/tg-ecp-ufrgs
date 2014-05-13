function Result = create_beatset(Database)
% Criaçao do conjunto de batidas

Records = fieldnames(Database);
for i = 1%1:numel(Records)
    ECG = Database.(Records{i});
    for j = 1%1:ECG.SignalCount
        Signal = ecgutilities.interpret(ECG,j);
        disp(['processing ' Signal.name '.' Signal.lead]);
        beatinfo = extract_beat_info(Signal, ECG.Annotations.atr);
        name = matlab.lang.makeValidName(Signal.name);
        lead = matlab.lang.makeValidName(Signal.lead);
        Result.(name).(lead) = beatinfo;
    end
end

function Result = extract_beat_info(Signal, atr)
import ecgutilities.*
[Beats,Rd,RR,Template] = ecgfilter.preprocess(Signal);
[idx1,idx2] = match_qrs(Signal.qrs, Rd, Signal.fs);
ST = extract_diagnosis(Signal.qrs(idx1), atr, 's', 'ST', Signal.id);
T = extract_diagnosis(Signal.qrs(idx1), atr, 'T', 'T', Signal.id);
Result.Fs = Signal.fs;
Result.R = Rd(idx2);
Result.RR = RR(idx2);
Result.Beats = Beats(:,idx2);
Result.Template = Template;
Result.ST = ST;
Result.T = T;