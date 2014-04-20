function Result = extract_rpeak_info(Signal,Fs,Rknown)

Filt = ecgfastcode.c_filter_double(Signal,Fs,3);
Rdetected = ecgfastcode.c_detect_qrs(Filt,Fs);
%[Rdetected,R2,TH1,TH2,RR] = ecgfastcode.c_detect_qrs(Filt,Fs);
%ecgutilities.plot_signal_qrs(Signal,Filt,Rdetected,R2,TH1,TH2,RR)
[A,B,R] = ecgutilities.merge_rpeaks(Rknown, Rdetected, Fs);
%ecgutilities.plot_signal_rcomp(Signal,A,B,R);
ecgmath.compute_statistics(A,B)

Result.A = A;
Result.B = B;
Result.R = R;