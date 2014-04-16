function Result = detect_qrs2(Signal, Fs)
%{
SignalI = ecgfilter.sogari_filter(Signal,Fs);
[R,R2,THs,THn,RRm] = ecgfeatures.sogari(SignalI,Fs);
figure;
hold on, grid on;
plot([SignalI THs THn]);
plot(R, SignalI(R), 'kx');
plot(R2, SignalI(R2), 'ko');
%figure, plot(RRm);
%}
SignalI = ecgfilter.sogari_filter(Signal,Fs);
Result = ecgfeatures.sogari_production(SignalI,Fs);