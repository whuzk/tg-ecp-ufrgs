% leitura do ecg
signal = utils.interpret_ecg(e0103,1);

% preprocessamento
tic
[R,RR,F,Beats,Templates,delay] = preprocess.preprocess(signal);
toc;

% extra�ao de caracteristicas do Rocha
tic;
[Rocha,RochaS1,RochaS2] = features.rocha_features(signal, RR, F, Beats);
toc;

% extra�ao de caracteristicas do Mohebbi
tic;
[Mohebbi,MohebbiS] = features.mohebbi_features(signal, F, Beats, Templates(:,30));
toc;

% extra�ao de caracteristicas do Gopalakrishnan
tic;
[Gopalak,GopalakS] = features.gopalak_features(RR, Beats);
toc;

% visualiza�ao
%plot.plot_rocha_features(Rocha, RochaS1, RochaS2);
%plot.plot_mohebbi_features(Mohebbi, MohebbiS);
%plot.plot_gopalak_features(Gopalak, GopalakS);