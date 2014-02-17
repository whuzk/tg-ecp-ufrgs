function net_train_rocha_networks

global Datasets RochaNet RochaFeatures RochaDiagnosis;

%{
%%
names = fieldnames(Datasets);
count = 0;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Rocha;
    count = count + size(Temp.F,1);
end

RochaFeatures = zeros(count,14);
RochaDiagnosis = zeros(count,2);
B = ones(1,10)/10;
endIndex = 0;
startIndex = 1;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Rocha;
    endIndex = endIndex + size(Temp.F,1);
    RochaFeatures(startIndex:endIndex,:) = filter(B,1,Temp.F);
    RochaDiagnosis(startIndex:endIndex,:) = 2*Temp.D-1;
    startIndex = endIndex + 1;
end
%}

%%
%{
a = RochaDiagnosis(:,1) > 0 | RochaDiagnosis(:,2) > 0;
indexNormalBeats = find(~a);
indexIschemicBeats = find(a);
nn = length(indexNormalBeats);
ni = length(indexIschemicBeats);
an = false(nn,1);
ai = false(ni,1);
an(randperm(nn,fix(0.7*nn))) = true;
ai(randperm(ni,fix(0.7*ni))) = true;
trainInd = [indexNormalBeats(an); indexIschemicBeats(ai)];
valInd = [indexNormalBeats(~an); indexIschemicBeats(~ai)];
trainInd = trainInd(randperm(end));
valInd = valInd(randperm(end));

net = feedforwardnet([13 6], 'trainlm');
net.layers{:}.transferFcn = 'tansig';
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd = valInd;
net.divideParam.testInd = [];
net.trainParam.epochs = 2000;
RochaNet = train(net, RochaFeatures', RochaDiagnosis');
%save('../resources/nnetworks.mat', 'RochaNet', '-append');
%}

%%
%{
%B = ones(1,10)/10;
%O = ecgrocha.ecg_classify_ischemic_beats(filter(B,1,Datasets.e0104.V4.Rocha.F));
%D = Datasets.e0104.V4.Rocha.D(:,1) | Datasets.e0104.V4.Rocha.D(:,2);
O = ecgrocha.ecg_classify_ischemic_beats(RochaFeatures);
D = RochaDiagnosis(:,1) > 0 | RochaDiagnosis(:,2) > 0;
Stats = utilities.compute_statistics(D, O);
disp(Stats);
%}
