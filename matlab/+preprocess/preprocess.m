function [R,RR,F,Beats,Template,delay] = preprocess(signal)

% obtem as amostras do sinal e a taxa de amostragem
data = signal.data - signal.inival;
Fs = signal.fs;

% filtragem do sinal
[sigI,sigD,sigM,sigN,delay] = mex.preprocess_filter(data, Fs);

% dete�ao de pontos e batidas
[R,RR,F,Beats,Template] = mex.preprocess_detect(sigI,sigD,sigM,sigN,Fs,30);