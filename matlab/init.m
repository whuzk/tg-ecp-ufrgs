% init.m
%   Inicializa?ao do ambiente de trabalho.
%
global methods leadnames edbleadnames
global stclasses tclasses classes measures

% nome dos metodos
methods = {
    'Rocha'
    'Mohebbi'
    'Gopalak'
};

% nomde das deriva?oes
leadnames = {
    'D1' 'D2' 'D3' 'MLI' 'MLII' 'MLIII' 'aVR' 'aVL' 'aVF' ...
    'V1' 'V2' 'V3' 'V4' 'V5' 'V6'
};

% nomde das derivacoes na base EDB
edbleadnames = leadnames([3:4 6 10:14]);

% nomde das classes
stclasses = {
    'Normal ST'
    'Elevated ST'
    'Depressed ST'
};
tclasses = {
    'Normal T'
    'Elevated T'
    'Depressed T'
};
classes = {
    [stclasses{1} ' and ' tclasses{1}]
    [stclasses{2} ' and ' tclasses{1}]
    [stclasses{3} ' and ' tclasses{1}]
    [stclasses{1} ' and ' tclasses{2}]
    [stclasses{2} ' and ' tclasses{2}]
    [stclasses{3} ' and ' tclasses{2}]
    [stclasses{1} ' and ' tclasses{3}]
    [stclasses{2} ' and ' tclasses{3}]
    [stclasses{3} ' and ' tclasses{3}]
};

% nome das estatisticas
measures = {
    'Sensitivity'
    'Specificity'
    'PositivePred'
    'NegativePred'
    'Accuracy'
    'FailureRate'
};