import ecgfastcode.*
%close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);
Signal = round(Signal*200);

Filt = c_filter_double(Signal,Fs);
[Result,R2,THs,THn,RRm] = c_detect_qrs(Filt,Fs);
%[Result2,R22,THs2,THn2,RRm2] = ecgfeatures.sogari_qrs(Filt,Fs);

%{
norm(Result-Result2)
norm(R2-R22)
norm(THs-THs2)
norm(THn-THn2)
norm(RRm-RRm2)
figure, plot(abs(Result-Result2))
figure, plot(abs(R2-R22))
figure, plot(abs(THs-THs2))
figure, plot(abs(THn-THn2))
figure, plot(abs(RRm-RRm2))
%}
%
figure;
hold on, grid on;
plot([Filt THs THn]);
plot(Result, Filt(Result), 'kx');
plot(R2, Filt(R2), 'ko');
figure, plot(RRm);
%