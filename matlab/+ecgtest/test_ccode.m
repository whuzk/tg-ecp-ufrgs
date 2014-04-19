import ecgfastcode.*
import ecgutilities.*
%close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(:,2);
Signal = round(Signal*200);

Filt = c_filter_double(Signal,Fs);
[R1,R2,TH1,TH2,RR] = c_detect_qrs(Filt,Fs);
%[R1,R2,TH1,TH2,RR] = ecgfeatures.sogari_qrs(Filt,Fs);
plot_signal_qrs(Filt,R1,R2,TH1,TH2,RR)