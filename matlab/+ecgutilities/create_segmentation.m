function Result = create_segmentation(Database)

Records = fieldnames(Database);
for i = 1:numel(Records)
    ECG = Database.(Records{i});
    for j = 1:ECG.SignalCount
        Signal = ecgutilities.interpret(ECG,j);
        disp(['processing ' Signal.name '.' Signal.lead]);
        qrsinfo = extract_qrs_info(Signal);
        name = matlab.lang.makeValidName(Signal.name);
        lead = matlab.lang.makeValidName(Signal.lead);
        Result.(name).(lead) = qrsinfo;
    end
end

function Result = extract_qrs_info(Signal)
data = Signal.data - Signal.inival;
[~,Rd,~,~,~,~,delay] = ecgfastcode.detect_qrs_double(data,Signal.fs);
%[Filt,Rd,R2,TH1,TH2,RR,delay] = ecgfastcode.detect_qrs_double(data,Signal.fs);
%ecgutilities.plot_signal_r(Signal,Rd-floor(delay));
%ecgutilities.plot_signal_qrs(Filt,Rd,R2,TH1,TH2,RR);
[A,B,R] = ecgutilities.merge_qrs(Signal.qrs, Rd-floor(delay), Signal.fs);
%ecgutilities.plot_signal_rcomp(Signal,A,B,R);
%ecgmath.compute_statistics(A,B)
Result.A = A;
Result.B = B;
Result.R = R;