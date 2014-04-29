import ecgutilities.*
import ecgfilter.*
import ecgmath.*
import ecgfastcode.*
%close all;

% load ecg
Signal = interpret(EDB.e0116,1);
data = Signal.data - Signal.inival;
dataInt = int16(data);
%{
% apply filters
Y1 = suppress_noise(Signal,Fs);
Y2 = baseline_filter(Signal,Fs);
Y3 = chu_filter(Signal,Fs);
Y4 = krishnan_filter(Signal,Fs);
Y5 = trahanias_filter(Signal,Fs);
Y6 = sun_filter(Signal,Fs);
[~,Y7] = tompkins_filter(Signal, Fs);

% plot figures
figure, plot(Signal);
figure, plot(Y1);
figure, plot(Y2);
figure, plot(Y3);
figure, plot(Y4);
figure, plot(Y5);
figure, plot(Y6);
figure, plot(Y7);
%}

tic;
[mmd1,lap1,delay1] = detect_fiducial_marks_double(data,Signal.lead,Signal.fs);
toc;
tic;
[mmd2,lap2,delay2] = detect_fiducial_marks_int(dataInt,Signal.lead,Signal.fs);
toc;
mmd2 = double(mmd2);
lap2 = double(lap2);
figure, plot([mmd1 mmd2]);
figure, plot([lap1 lap2]);
rms(mmd1-mmd2)
rms(lap1-lap2)

%data = data(1:50000);
%FrameSize = 2*floor(Signal.fs*0.6)+1;
%[Beats,R,RR] = extract_beats(data,Signal.lead,Signal.fs,FrameSize);
%figure, plot(Beats(:,1));
%figure, plot(Beats(:,100));
%figure, plot(Beats(:,1000));