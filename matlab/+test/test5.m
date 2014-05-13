% leitura do ecg
signal = utils.interpret(e0103,1);

% preprocessamento
tic
[R,RR,F,Beats,Templates,delay] = preprocess.preprocess(signal);
toc;

% extraçao de caracteristicas do Rocha
tic;
[Rocha,RochaS1,RochaS2] = features.rocha_features(signal, RR, F, Beats);
toc;

% extraçao de caracteristicas do Mohebbi
tic;
[Mohebbi,MohebbiS] = features.mohebbi_features(signal, F, Beats, Templates(:,30));
toc;

% extraçao de caracteristicas do Gopalakrishnan
tic;
[Gopalak,GopalakS] = features.gopalak_features(RR, Beats);
toc;

% visualizaçao
%utils.plot_rocha(Rocha, RochaS1, RochaS2);
%utils.plot_mohebbi(Mohebbi, MohebbiS);
%utils.plot_gopalak(Gopalak, GopalakS);