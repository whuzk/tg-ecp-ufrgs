% init.m
%   Inicializa�ao do ambiente de trabalho.
%
global leadnames measures methods

% nomde das deriva�oes
leadnames = {
    'D1' 'D2' 'D3' 'MLI' 'MLII' 'MLIII' 'aVR' 'aVL' 'aVF' ...
    'V1' 'V2' 'V3' 'V4' 'V5' 'V6'
};

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

% registros de ecg
EDB = load('C:\physiobank\database\edb.mat', 'e0119');
%MITDB = load('C:\physiobank\database\mitdb.mat');
