% init.m
%   Inicializaçao do ambiente de trabalho.
%

% nome da derivaçao
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
EDB = load('C:\physiobank\database\edb.mat', 'e0116');
%MITDB = load('C:\physiobank\database\mitdb.mat');
%QTDB = load('C:\physiobank\database\qtdb.mat');

% inicializa gerador pseudo-randomico
%rng('shuffle');
