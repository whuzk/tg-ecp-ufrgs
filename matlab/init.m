% init.m
%   Inicializaçao do ambiente de trabalho.
%
global leadnames measures

% nomde das derivaçoes
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

% registros de ecg
x100 = load('C:\physiobank\database\mitdb\100.mat');
e0103 = load('C:\physiobank\database\edb\e0103.mat');
e0104 = load('C:\physiobank\database\edb\e0104.mat');
e0105 = load('C:\physiobank\database\edb\e0105.mat');