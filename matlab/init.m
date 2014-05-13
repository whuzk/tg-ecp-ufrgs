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
x100 = load('matfiles\x100.mat');
x105 = load('matfiles\x105.mat');
e0103 = load('matfiles\e0103.mat');
e0116 = load('matfiles\e0116.mat');