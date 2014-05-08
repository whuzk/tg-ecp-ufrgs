import ecgutilities.*

% le uma das derivaçoes do ecg
signal = interpret(EDB.e0119,2);
data = signal.data - signal.inival;

% detecta os picos de onda R
[sigI,delay] = mex.tompkins_filter(data, signal.fs);
[qrs1,qrs2,th1,th2,rr1,rr2] = mex.detect_qrs(sigI, signal.fs);

% compara com as anotaçoes
Radj = adjust_qrs(data,signal.lead,qrs1-floor(delay),signal.fs);
[a,b,c] = merge_qrs(signal.qrs, Radj, signal.fs);
stats = ecgmath.compute_statistics(a,b);

% plota os graficos
plot_rpeaks(data, Radj);
plot_thresholds(sigI,qrs1,qrs2,th1,th2,rr1,rr2);
plot_comparison(data,a,b,c);

% mostra o resultado da comparaçao
ds = dataset([{stats'} measures{:}]);
ds.RecordName = signal.name;
disp(ds);