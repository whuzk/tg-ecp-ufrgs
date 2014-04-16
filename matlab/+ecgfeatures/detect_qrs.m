function Result = detect_qrs(Signal, Fs)
%{
[SignalF,SignalI] = ecgfilter.tompkins_filter(Signal, Fs);
[R,R2,THs,THn,THs2,THn2] = ...
    ecgfeatures.tompkins_adapted(SignalF, SignalI, Fs);
figure;
hold on, grid on;
plot([SignalF THs2 THn2]);
figure;
hold on, grid on;
plot([SignalI THs THn]);
plot(R, SignalI(R), 'kx');
plot(R2, SignalI(R2), 'ko');
%}
[SignalF,SignalI] = ecgfilter.tompkins_filter(Signal, Fs);
Result = ecgfeatures.tompkins_production(SignalF, SignalI, Fs);