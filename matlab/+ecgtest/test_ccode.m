import ecgfilter.*
%close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);
Signal = round(Signal*200);
SignalInt32 = int32(Signal);

% filter
%{
h = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/4;
Y1 = wconv1(Signal, h, 'same');
Y2 = ecgfilter.mobd(Y1,3);
Ws = floor(0.05*Fs)*2+1;
Y3 = wconv1(Y2,ones(1,Ws)/Ws,'same');
%}

Z3 = double(c_filter(SignalInt32,Fs));
Z3 = c_filter_double(Signal,Fs);

%figure, plot(Signal);
norm(Y3-Z3)
figure, plot([Y3 Z3]);