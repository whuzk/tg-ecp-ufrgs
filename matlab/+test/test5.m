% leitura do ecg
record = e0104;
recordIndex = 1;
signal = utils.interpret_ecg(record, recordIndex);

% preprocessamento
tic
[R,RR,F,Beats,Template,delay] = preprocess.preprocess(signal);
toc;

% extraçao de caracteristicas do Rocha
tic;
[Rocha,RochaS1,RochaS2] = features.rocha_features(signal, RR, F, Beats);
toc;

% extraçao de caracteristicas do Mohebbi
tic;
[Mohebbi,MohebbiS] = features.mohebbi_features(signal, F, Beats, Template);
toc;

% extraçao de caracteristicas do Gopalakrishnan
tic;
[Gopalak,GopalakS] = features.gopalak_features(signal, RR, Beats);
toc;

% visualizaçao das caracteristicas
%plot.plot_rocha_features(Rocha, RochaS1, RochaS2);
%plot.plot_mohebbi_features(Mohebbi, MohebbiS);
%plot.plot_gopalak_features(Gopalak, GopalakS);

% extraçao do diagnostico de isquemia
[STdiag,Tdiag,idx] = utils.extract_ischemia_diag(record, signal.fs, R-ceil(delay), recordIndex-1);

% composiçao do conjunto de dados para treinamento das redes neurais
Datasets.Rocha = [Rocha(:,idx)' STdiag Tdiag];
Datasets.Mohebbi = [Mohebbi(:,idx)' STdiag Tdiag];
Datasets.Gopalak = [Gopalak(:,idx)' STdiag Tdiag];

% composiçao das informaçoes basicas extraidas no preprocessamento
Preprocess.R = R(idx)-ceil(delay);
Preprocess.RR = RR(idx);
Preprocess.F = F(:,idx);
Preprocess.Beats = Beats(:,idx);
Preprocess.Template = Template;
Preprocess.delay = delay;

% salvamento em arquivo
filename = ['extracted_' signal.name '.mat'];
save(filename, 'Preprocess', 'Datasets');
