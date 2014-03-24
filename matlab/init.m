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
};

% nome dos registros
load('.\resources\dbnames.mat', 'EDBnames');

% registros de ecg
Database = load('C:\ecg\edb.mat');
%Database = load('C:\ecg\edb.mat', EDBnames{1:10});

% inicializa gerador pseudo-randomico
rng('shuffle');
