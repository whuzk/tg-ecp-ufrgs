% init.m
%   Inicializaçao do ambiente de trabalho.
%
global leadnames measures methods
global pwavepolarity rwavepolarity twavepolarity

% nomde das derivaçoes
leadnames = {
    'D1' 'D2' 'D3' 'MLI' 'MLII' 'MLIII' 'aVR' 'aVL' 'aVF' ...
    'V1' 'V2' 'V3' 'V4' 'V5' 'V6'
};

% polaridade das ondas
pwavepolarity = containers.Map(leadnames, ...
    [1  1  1  1  1  1 -1  1  1  1  1  1  1  1  1]);
rwavepolarity = containers.Map(leadnames, ...
    [1  1 -1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1]);
twavepolarity = containers.Map(leadnames, ...
    [1  1  1  1  1 -1 -1 -1  1  1  1  1  1  1  1]);

% nome das estatisticas
measures = {
    'Sensitivity'
    'Specificity'
    'PositivePred'
    'NegativePred'
    'Accuracy'
    'FalseDetection'
};

% nome dos metodos
methods = {'Rocha', 'Mohebbi', 'Gopalakrishnan'};

% nome dos registros
load('.\matfiles\dbnames.mat');

% registros de ecg
EDB = load('C:\physiobank\database\edb.mat');
%MITDB = load('C:\physiobank\database\mitdb.mat');
