% init.m
%   Inicializašao do ambiente de trabalho.
%

% nome da derivašao
LeadName = 'V4';

% nome dos metodos
Methods = {'Rocha', 'Mohebbi', 'Gopalakrishnan'};

% nome das estatisticas
Measures = {
    'Sensitivity'
    'Specificity'
    'PositivePred'
    'NegativePred'
    'Accuracy'
    'False Detection Rate'
};

% nome dos registros
load('.\matfiles\dbnames.mat');

% registros de ecg
EDB = load('C:\ecg\edb.mat', 'e0116');
%MITDB = load('C:\ecg\mitdb.mat');
%QTDB = load('C:\ecg\qtdb.mat');

% inicializa gerador pseudo-randomico
%rng('shuffle');
