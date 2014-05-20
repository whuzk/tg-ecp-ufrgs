function [STdiag,Tdiag,idx] = extract_ischemia_diag(ecg, Fs, R, sid)

atr = ecg.Annotations.atr;
[A,idx] = utils.match_qrs(atr.Sample, R, Fs);
R = atr.Sample(A);

ST = utils.extract_diagnosis(R, atr, 's', 'ST', sid);
T = utils.extract_diagnosis(R, atr, 'T', 'T', sid);

STdiag = ST.diagVal;
Tdiag = T.diagVal;