function [Info,leadname] = extract_datasets(record, recordIndex)

% leitura do ecg
signal = utils.interpret_ecg(record, recordIndex);
leadname = signal.lead;

% preprocessamento
[R,RR,F,Beats,Template,delay] = preprocess.preprocess(signal);

% extra�ao de caracteristicas do Rocha
Rocha = features.rocha_features(signal, RR, F, Beats);

% extra�ao de caracteristicas do Mohebbi
Mohebbi = features.mohebbi_features(signal, F, Beats, Template);

% extra�ao de caracteristicas do Gopalakrishnan
Gopalak = features.gopalak_features(signal, RR, Beats);

% extra�ao do diagnostico de isquemia
[STdiag,Tdiag,idx] = utils.extract_ischemia_diag(record, signal.fs, R-ceil(delay), recordIndex-1);

% composi�ao do conjunto de dados para treinamento das redes neurais
Info.Datasets.Rocha = [Rocha(:,idx)' STdiag Tdiag];
Info.Datasets.Mohebbi = [Mohebbi(:,idx)' STdiag Tdiag];
Info.Datasets.Gopalak = [Gopalak(:,idx)' STdiag Tdiag];

% composi�ao das informa�oes basicas extraidas no preprocessamento
Info.Preprocess.R = R(idx)-ceil(delay);
Info.Preprocess.RR = RR(idx);
Info.Preprocess.F = F(:,idx);
Info.Preprocess.Beats = Beats(:,idx);
Info.Preprocess.Template = Template;
Info.Preprocess.delay = delay;