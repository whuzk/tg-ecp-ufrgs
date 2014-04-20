function Result = extract_beat_info(Signal,Fs,Rknown,id,Annot)

[Beats,R,RR,Template] = ecgfilter.preprocess(Signal, Fs);
[index1,index2] = ecgutilities.find_close_beats(Rknown, R, Fs);
STseg = extract_diagnosis(Annot, Bp(index1), 's', 'ST', id);
Twave = extract_diagnosis(Annot, Bp(index1), 'T', 'T', id);

Result.Fs = Fs;
Result.R = R(index2);
Result.RR = RR(index2);
Result.Beats = Beats(:,index2);
Result.Template = Template;
Result.STseg = STseg;
Result.Twave = Twave;