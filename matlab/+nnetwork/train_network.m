function [net,tr] = train_network(inputs, targets, layers, trainfcn, transfcn)

net = feedforwardnet(layers, trainfcn);
net.layers{:}.transferFcn = transfcn;

net.divideFcn = 'dividerand';
net.divideMode = 'sample';
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.15;

net.trainParam.showWindow = false;
net.trainParam.showCommandLine = true;
net.trainParam.show = 20;
net.trainParam.time = 5*60;
net.trainParam.epochs = 1000;
net.trainParam.max_fail = 10;
net.trainParam.min_grad = 1e-10;

net.plotFcns = {'plotconfusion', 'plotperform', 'plottrainstate', ...
    'ploterrhist', 'plotregression'};

[net,tr] = train(net,inputs,targets);