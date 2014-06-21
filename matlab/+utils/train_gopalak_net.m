basedir = '/Users/diegosogari/Documents/MATLAB/';
dataset = load([basedir 'gopalak_dataset.mat']);
x = dataset.V1(:,1:end-2)';
t = dataset.V1(:,end-1:end)';
t = 2*t-1;

%%
net = feedforwardnet(20);
net.divideFcn = 'dividerand';
net.divideMode = 'sample';
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.15;
net.trainParam.epochs = 2000;
net.trainFcn = 'trainlm';
net.trainParam.min_grad = 1e-10;
net.trainParam.max_fail = 20;
net.plotFcns = {'plotconfusion', 'plotperform', 'plottrainstate', ...
    'ploterrhist', 'plotregression', 'plotfit'};
[net,tr] = train(net,x,t);

%%
y = net(x);
%e = gsubtract(t,y);
%tind = vec2ind(t);
%yind = vec2ind(y);
%percentErrors = sum(tind ~= yind)/numel(tind);
performance = perform(net,t,y)

%%
trainTargets = t .* tr.trainMask{1};
valTargets = t  .* tr.valMask{1};
testTargets = t  .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y)
valPerformance = perform(net,valTargets,y)
testPerformance = perform(net,testTargets,y)

%%
Fs = 250;
gensim(net,Fs);