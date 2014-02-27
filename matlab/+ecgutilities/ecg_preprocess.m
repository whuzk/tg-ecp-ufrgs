function [Beats,Rpeaks,RR,Template] = ecg_preprocess(Signal, Fs)

MaxBeatLength = 2*fix(Fs*0.6)+1;

% remoçao de ruido
Signal = suppress_noise(Signal,Fs);

% segmentaçao e extraçao das batidas
Rpeaks = ecgutilities.ecg_segment(Signal,Fs);
RR = diff(Rpeaks(1:end-1));
Rpeaks = Rpeaks(2:end-1);
L = min(RR,MaxBeatLength);
Beats = extract_beats(Signal,Rpeaks,L);
%ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;

% remoçao da linha de base e enquadramento
tic;
m = length(Rpeaks);
for i = 1:m
    Beat = remove_baseline(Beats{i}, Fs);
    Beats{i} = frame_beat(Beat,MaxBeatLength);
end
toc;

% construçao de template
Tc = 30;
FirstBeats = cell2mat(Beats(1:Tc));
Template = mean(FirstBeats,2);
%figure, plot(Template);
%pause;

% remoçao de batidas anomalas
index = false(1,m);
for i = Tc+1:m
    Correlation = corr(Beats{i}, Template);
    %if Correlation < 0.75
    %    figure(10), plot(Beats{i});
    %    Correlation
    %    pause;
    %end
    index(i) = Correlation >= 0.75;
end
Beats = cell2mat(Beats(index));
Rpeaks = Rpeaks(index);
%ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;

function Result = suppress_noise(Signal, Fs)
Wn = 40 * 2/Fs;
[B,A] = butter(4, Wn);
Result = filter(B, A, Signal);

function Result = extract_beats(Signal, Rpeaks, RR)
m = length(Rpeaks);
half = fix(RR/2);
Result = cell(1,m);
for i = 1:m
    left = Rpeaks(i)-half(i);
    right = Rpeaks(i)+half(i);
    Result{i} = Signal(left:right);
end

function Result = remove_baseline(Signal, Fs)
n = length(Signal);
Rpeak = fix(n/2)+1;
[x,y] = ecgutilities.ecg_baseline_points(Signal, Fs, Rpeak);
Result = Signal - spline(x,y,(1:n)');

function Result = frame_beat(Signal, l)
D = l-length(Signal);
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    Result = Signal(1+a:end-b);
else
    Result = [zeros(a,1); Signal; zeros(b,1)];
end
