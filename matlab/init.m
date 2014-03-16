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

% registros de ecg
Database = load('C:\ecg\edb.mat', ...
    'e0103', ...
    'e0104', ...
    'e0105', ...
    'e0107', ...
    'e0113', ...
    'e0119', ...
    'e0147' ...
);

% inicializa gerador pseudo-randomico
rng('shuffle');
