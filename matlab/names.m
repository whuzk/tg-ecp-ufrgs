% these variables contain names of the ECG records for which the
% segmentation algorithm by Pan-Tompkins falls short of 99.0% in
% terms of sensibility and/or positive predictivity
D3 = {};
MLI = {};
MLIII = {
    'e0111' 'e0112' 'e0116' 'e0118' ...
    'e0119' 'e0121' 'e0129' 'e0139' ...
    'e0151' 'e0155' 'e0159' 'e0605' ...
    'e0611' 'e0613'
};
V1 = {'e0614'};
V2 = {'e0305' 'e0413' 'e0415'};
V3 = {};
V4 = {'e0116'};
V5 = {'e0413' 'e0614' 'e0801' 'e0817'};