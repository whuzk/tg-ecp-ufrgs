function Result = extract_rpeak_info(Signal,Fs,Rknown)

[Filt,delay] = ecgfastcode.filter_double(Signal,Fs,3);

Rdetected = ecgfastcode.detect_qrs(Filt,Fs);
%[Rdetected,R2,TH1,TH2,RR] = ecgfastcode.detect_qrs(Filt,Fs);
%ecgutilities.plot_signal_r(Signal,Rdetected-delay);
%ecgutilities.plot_signal_qrs(Filt,Rdetected,R2,TH1,TH2,RR);

[A,B,R] = ecgutilities.merge_rpeaks(Rknown, Rdetected-delay, Fs);
%ecgutilities.plot_signal_rcomp(Signal,A,B,R);
ecgmath.compute_statistics(A,B)

Result.A = A;
Result.B = B;
Result.R = R;