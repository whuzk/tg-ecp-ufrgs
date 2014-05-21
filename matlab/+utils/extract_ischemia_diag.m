function [STdiag,Tdiag,idx2] = extract_ischemia_diag(ecg, Fs, R, sid)

atr = ecg.Annotations.atr;
[idx1,idx2] = utils.match_qrs(atr.Sample, R, Fs);
Rannot = atr.Sample(idx1);

ST = utils.extract_diagnosis(Rannot, atr, 's', 'ST', sid);
T = utils.extract_diagnosis(Rannot, atr, 'T', 'T', sid);

STdiag = ST.diagVal;
Tdiag = T.diagVal;