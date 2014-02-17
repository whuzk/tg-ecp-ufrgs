function Result = train_network(MethodDataset)
% Treinamento das redes neurais de um metodo

%
B = ones(1,10)/10;
F = filter(B,1,MethodDataset.Features);
D = 2*MethodDataset.Diagnosis-1;

%
a = D(:,1) > 0 | D(:,2) > 0;
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

%
net = feedforwardnet(20, 'trainlm');
net.layers{:}.transferFcn = 'tansig';
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd = valInd;
net.divideParam.testInd = [];
net.trainParam.epochs = 2000;
Result = train(net, F', D');
