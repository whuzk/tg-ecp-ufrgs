% init.m
%   Inicializaçao do ambiente de trabalho.
%

% nome da derivaçao
LeadName = 'V4';
%LeadName = 'MLIII';

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
Database = load('C:\ecg\edb.mat', EDBnames{1:10});

% inicializa gerador pseudo-randomico
rng('shuffle');
