% these variables contain names of the ECG records for which the
% segmentation algorithm by Pan-Tompkins falls short of 99.0% in
% terms of sensibility and/or positive predictivity

%% EDB
EDB.D3 = {}; %none
EDB.MLI = {}; %none
EDB.MLIII = {
    'e0108' 'e0111' 'e0112' 'e0116' ...
    'e0118' 'e0119' 'e0129' 'e0151' ...
    'e0155' 'e0159' 'e0611'
};
EDB.V1 = {}; %none
EDB.V2 = {'e0305' 'e0415'};
EDB.V3 = {}; %none
EDB.V4 = {'e0116'};
EDB.V5 = {'e0801' 'e0817'};

%% MITDB
MITDB.MLII = {
    'e105' 'e106' 'e108' 'e201' 'e203' ...
    'e207' 'e210' 'e219' 'e231' 'e232'
};
MITDB.V1 = {
    'e101' 'e105' 'e106' 'e108' 'e115' ...
    'e200' 'e201' 'e203' 'e207' 'e208' ...
    'e209' 'e219' 'e231' 'e232'
};
MITDB.V2 = {'e103' 'e104'};
MITDB.V4 = {}; %none
MITDB.V5 = {}; %none

%% QTDB
QTDB.MLII = {}; %none
QTDB.V1 = {}; %none
QTDB.V2 = {}; %none
QTDB.V3 = {'sele0129'};
QTDB.V5 = {}; %none
QTDB.D3 = {}; %none
QTDB.D4 = {}; %none
QTDB.CM2 = {}; %none
QTDB.CM5 = {}; %none
QTDB.CC5 = {}; %none
QTDB.ML5 = {}; %none
QTDB.V1V2 = {}; %none
QTDB.V2V3 = {}; %none
QTDB.V4V5 = {}; %none
