% leitura do ecg
record = load('C:\physiobank\database\edb\e0103.mat');
signal = utils.interpret_ecg(record, 1);

% preprocessamento
tic
[R,RR,F,Beats,Template,delay] = preprocess.preprocess(signal);
toc;

% extrašao de caracteristicas do Rocha
tic;
[Rocha,RochaS1,RochaS2] = features.rocha_features(signal, RR, F, Beats);
toc;

% extrašao de caracteristicas do Mohebbi
tic;
[Mohebbi,MohebbiS] = features.mohebbi_features(signal, F, Beats, Template);
toc;

% extrašao de caracteristicas do Gopalakrishnan
tic;
[Gopalak,GopalakS] = features.gopalak_features(signal, RR, Beats);
toc;

% visualizašao das caracteristicas
%plot.plot_rocha_features(Rocha, RochaS1, RochaS2);
%plot.plot_mohebbi_features(Mohebbi, MohebbiS);
%plot.plot_gopalak_features(Gopalak, GopalakS);