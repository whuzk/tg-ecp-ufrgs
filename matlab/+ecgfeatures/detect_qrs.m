function Result = detect_qrs(Signal, Fs)

Filt = ecgfastcode.c_filter_double(Signal,Fs);
[Result,R2,THs,THn,RRm] = ecgfeatures.sogari_qrs(Filt,Fs);

figure;
hold on, grid on;
plot([Filt THs THn]);
plot(Result, Filt(Result), 'kx');
plot(R2, Filt(R2), 'ko');
figure, plot(RRm);