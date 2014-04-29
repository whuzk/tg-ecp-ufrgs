import ecgfastcode.*
import ecgutilities.*
import ecgfilter.*
%close all;
%{
% load ecg
Signal = interpret(EDB.e0116,2);
data = Signal.data - Signal.inival;
%
% detect qrs
[R,RR,delay] = prod_detect_qrs_double(data,Signal.fs,50,3);
[data,delay2] = suppress_noise(data,Signal.fs);
R = R - floor(delay - delay2);
R = adjust_qrs(data,Signal.lead,R,Signal.fs);
plot_signal_r(data, R);
%}
% extract beats
FrameSize = 2*floor(Signal.fs*0.6)+1;
Beats = extract_beats(data,Signal.lead,Signal.fs,R,FrameSize);
%{
plot_signal(Beats(:,1), 'Beat #1');
plot_signal(Beats(:,100), 'Beat #100');
plot_signal(Beats(:,1000), 'Beat #1000');
%}