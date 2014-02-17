function net_train_mohebbi_network
% Treinamento da rede neural do metodo do Mohebbi

global Datasets MohebbiNet MohebbiFeatures MohebbiDiagnosis;

%{
%%
names = fieldnames(Datasets);
count = 0;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Mohebbi;
    count = count + size(Temp.F,1);
end

MohebbiFeatures = zeros(count,20);
MohebbiDiagnosis = zeros(count,2);
B = ones(1,10)/10;
endIndex = 0;
startIndex = 1;
for i = 1:numel(names)
    Temp = Datasets.(names{i}).V4.Mohebbi;
    endIndex = endIndex + size(Temp.F,1);
    MohebbiFeatures(startIndex:endIndex,:) = filter(B,1,Temp.F);
    MohebbiDiagnosis(startIndex:endIndex,1) = 2*Temp.D-1;
    MohebbiDiagnosis(startIndex:endIndex,2) = 2*(~Temp.D)-1;
    startIndex = endIndex + 1;
end
%}

%{
%%
indexNormalBeats = find(MohebbiDiagnosis(:,2) > 0);
indexIschemicBeats = find(MohebbiDiagnosis(:,1) > 0);
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

%
net = feedforwardnet(20, 'trainlm');
net.layers{:}.transferFcn = 'tansig';
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd = valInd;
net.divideParam.testInd = [];
net.trainParam.epochs = 2000;
MohebbiNet = train(net, MohebbiFeatures', MohebbiDiagnosis');
%save('../resources/nnetworks.mat', 'MohebbiNet', '-append');
%}

%{
%%
%B = ones(1,10)/10;
%O = ecgmohebbi.ecg_classify_ischemic_beats(filter(B,1,Datasets.e0104.V4.Mohebbi.F));
%D = Datasets.e0104.V4.Mohebbi.D;
O = ecgmohebbi.ecg_classify_ischemic_beats(MohebbiFeatures);
D = MohebbiDiagnosis(:,1) > 0;
Stats = utilities.compute_statistics(D, O);
disp(Stats);
%}