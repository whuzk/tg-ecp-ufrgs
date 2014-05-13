% le uma das derivaçoes do ecg
signal = utils.interpret(x100,1);
data = signal.data - signal.inival;

% detecta os picos de onda R usando codigo c
[sigI,delay] = mex.tompkins_filter(data, signal.fs);
[qrs1,qrs2,th1,th2,rr1,rr2] = mex.detect_qrs(sigI, signal.fs);

% compara com as anotaçoes
Radj = qrs1 - floor(delay);
[a,b,c] = utils.merge_qrs(signal.qrs, Radj, signal.fs);
stats = utils.compute_statistics(a,b);

% plota os graficos
utils.plot_rpeaks(data, Radj);
utils.plot_thresholds(sigI,qrs1,qrs2,th1,th2,rr1,rr2);
utils.plot_comparison(data,a,b,c);

% mostra o resultado da comparaçao
ds = dataset([{stats'} measures{:}]);
ds.RecordName = signal.name;
disp(ds);