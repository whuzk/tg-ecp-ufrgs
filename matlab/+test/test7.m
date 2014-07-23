record = load('C:\physiobank\database\edb\e0161.mat');
signal = utils.interpret_ecg(record, 2);
data = signal.data - signal.inival;
Fs = signal.fs;

[sigI,sigD,sigM,sigN,delay] = mex.preprocess_filter(data, Fs);
[R,RR,F,Beats,Template] = mex.preprocess_detect(sigI,sigD,sigM,sigN,Fs,30);
%R = mex.detect_qrs(sigI, Fs);
%F2 = preprocess.fiducial_marks(sigD,sigM,R,Fs);
%plot.plot_fiducial_marks(sigD,F2);
plot.plot_fiducial_marks(sigM,F2);
plot.plot_fiducial_marks(sigN,F2);