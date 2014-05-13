function net_train_gopalak_network
% Treinamento da rede neural do metodo do Mohebbi

global Datasets GopalakNets GopalakFeatures GopalakDiagnosis;

%{
%%
names = fieldnames(Datasets);
count = 0;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Gopalak;
    count = count + size(Temp.F,1);
end

GopalakFeatures = zeros(count,50);
GopalakDiagnosis = zeros(count,2);
B = ones(1,10)/10;
endIndex = 0;
startIndex = 1;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Gopalak;
    endIndex = endIndex + size(Temp.F,1);
    GopalakFeatures(startIndex:endIndex,:) = filter(B,1,Temp.F);
    GopalakDiagnosis(startIndex:endIndex,:) = 2*Temp.D-1;
    startIndex = endIndex + 1;
end
%}

%{
%%
a = GopalakDiagnosis(:,1) > 0 | GopalakDiagnosis(:,2) > 0;
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

net = feedforwardnet(20, 'trainlm');
net.layers{:}.transferFcn = 'tansig';
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd = valInd;
net.divideParam.testInd = [];
net.trainParam.epochs = 2000;
GopalakNets{5} = train(net, GopalakFeatures', GopalakDiagnosis');
%save('../resources/nnetworks.mat', 'GopalakNets', '-append');
%}

%%
%{
%B = ones(1,10)/10;
%O = ecggopalak.ecg_classify_ischemic_beats(filter(B,1,Datasets.e0104.V4.Gopalak.F));
%D = Datasets.e0104.V4.Gopalak.D(:,1) | Datasets.e0104.V4.Gopalak.D(:,2);
O = ecggopalak.ecg_classify_ischemic_beats(GopalakFeatures);
D = GopalakDiagnosis(:,1) > 0 | GopalakDiagnosis(:,2) > 0;
Stats = utilities.compute_statistics(D, O);
disp(Stats);
%}