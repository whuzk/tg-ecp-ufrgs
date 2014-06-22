function Result = train_network(dataset, layers, trainfcn)

x = dataset(:,1:end-2)';
t1 = dataset(:,end-1)' ~= 0;
t2 = dataset(:,end)' ~= 0;

net = feedforwardnet(layers);
%net.layers{:}.transferFcn = 'tansig';
net.divideFcn = 'dividerand';
net.divideMode = 'sample';
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.15;
net.trainFcn = trainfcn;
net.trainParam.showWindow = false;
net.trainParam.showCommandLine = true;
net.trainParam.show = 20;
net.trainParam.epochs = 100;
net.trainParam.min_grad = 1e-10;
net.trainParam.max_fail = 10;
net.plotFcns = {'plotconfusion', 'plotperform', 'plottrainstate', ...
    'ploterrhist', 'plotregression'};

[net1,tr1] = train(net,x,t1);
[net2,tr2] = train(net,x,t2);
Result.nets = [net1 net2];
Result.recs = [tr1 tr2];