% le uma das deriva�oes do ecg
recordIndex = 1;
signal = utils.interpret_ecg(e0104, recordIndex);
data = signal.data - signal.inival;
Fs = signal.fs;

% filtragem
[sigI,sigD,sigM,sigN,delay] = mex.preprocess_filter(data, Fs);
[qrs1,qrs2,th1,th2,rr1,rr2] = mex.detect_qrs(sigI, Fs);

% compara com as anota�oes
Radj = qrs1 - floor(delay);
[a,b,c] = utils.merge_qrs(signal.qrs, Radj, Fs);
stats = utils.compute_statistics(a,b);

% plota os graficos
plot.plot_rpeaks(data, Radj);
plot.plot_qrs_detection(sigI,qrs1,qrs2,th1,th2,rr1,rr2);
plot.plot_rpeak_comparison(data,a,b,c);

% mostra o resultado da compara�ao
ds = dataset([{stats'} measures{:}]);
ds.RecordName = signal.name;
disp(ds);