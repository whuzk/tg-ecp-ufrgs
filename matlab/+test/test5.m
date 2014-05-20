% leitura do ecg
recordIndex = 1;
signal = utils.interpret_ecg(e0103, recordIndex);

% preprocessamento
tic
[R,RR,F,Beats,Template,delay] = preprocess.preprocess(signal);
toc;

% extra�ao de caracteristicas do Rocha
tic;
[Rocha,RochaS1,RochaS2] = features.rocha_features(signal, RR, F, Beats);
toc;

% extra�ao de caracteristicas do Mohebbi
tic;
[Mohebbi,MohebbiS] = features.mohebbi_features(signal, F, Beats, Template);
toc;

% extra�ao de caracteristicas do Gopalakrishnan
tic;
[Gopalak,GopalakS] = features.gopalak_features(signal, RR, Beats);
toc;

% extra�ao do diagnostico de isquemia
[STdiag,Tdiag,idx] = utils.extract_ischemia_diag(e0103, signal.fs, R-ceil(delay), recordIndex-1);

% visualiza�ao das caracteristicas
%plot.plot_rocha_features(Rocha, RochaS1, RochaS2);
%plot.plot_mohebbi_features(Mohebbi, MohebbiS);
%plot.plot_gopalak_features(Gopalak, GopalakS);

% composi�ao do conjunto de dados para treinamento das redes neurais
Datasets.Rocha = [Rocha' STdiag Tdiag];
Datasets.Mohebbi = [Mohebbi' STdiag Tdiag];
Datasets.Gopalak = [Gopalak' STdiag Tdiag];

% composi�ao das informa�oes basicas extraidas no preprocessamento
Preprocess.R = R(idx);
Preprocess.RR = RR(idx);
Preprocess.F = F(:,idx);
Preprocess.Beats = Beats(:,idx);
Preprocess.Template = Template;
Preprocess.delay = delay;

% salvamento em arquivo
filename = ['extracted_' signal.name '.mat'];
save(filename, 'Preprocess', 'Datasets');
