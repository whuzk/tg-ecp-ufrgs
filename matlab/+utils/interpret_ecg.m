function Result = interpret_ecg(Record,RecordIndex)

idx = ismember(Record.Annotations.atr.Type, 'NLRBAaJSVFejnE/fQ?');

Result.id = Record.Info(RecordIndex).SignalIndex;
Result.fs = sscanf(Record.Info(RecordIndex).SamplingFrequency,'%d');
Result.name = sscanf(Record.Info(RecordIndex).RecordName,'%[^.]');
Result.lead = Record.Info(RecordIndex).Description;
Result.res = sscanf(Record.Info(RecordIndex).AdcResolution,'%d');
Result.gain = sscanf(Record.Info(RecordIndex).Gain,'%d');
Result.inival = Record.Info(RecordIndex).InitialValue;
Result.zero = Record.Info(RecordIndex).AdcZero;
Result.len = Record.Info(RecordIndex).LengthSamples;
Result.time = (0:Result.len-1)'./Result.fs;
Result.data = Record.SignalData(:,RecordIndex);
Result.qrs = Record.Annotations.atr.Sample(idx);