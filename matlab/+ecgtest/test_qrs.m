dbpath = 'I:\AppData\physiobank\database\edb';
record = 'e0103';
old_path = pwd;
cd(dbpath);

gqrs(record,[],[],1,[]);
%sqrs(record,[],[],1,[]);
%wqrs(record,[],[],1,[]);
ann = rdann(record,'qrs',1,[],[]);
[tm,sig] = rdsamp(record,1,[],[]);

cd(old_path);

ecgutilities.plot_signal_r(Signal,ann);
[A,B,R] = ecgutilities.merge_rpeaks(Bp, ann, Fs);
ecgutilities.plot_signal_rcomp(Signal,A,B,R);
ecgmath.compute_statistics(A,B)
