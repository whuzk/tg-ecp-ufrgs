function [Beats,Rpeaks,RR,Template] = ecg_preprocess(Signal, Fs)

% remoçao de ruido
Signal = suppress_noise(Signal,Fs);

% segmentaçao e extraçao das batidas
Rpeaks = ecgutilities.ecg_segment(Signal,Fs);
RR = diff(Rpeaks(1:end-1));
Rpeaks = Rpeaks(2:end-1);
Beats = extract_beats(Signal,Rpeaks,RR);

% remoçao da linha de base e enquadramento
m = length(Rpeaks);
n = 2*fix(Fs*0.6)+1;
for i = 1:m
    Beat = remove_baseline(Beats{i}, Fs);
    Beats{i} = frame_beat(Beat,n);
end

% construçao de template
Tc = 30;
FirstBeats = cell2mat(Beats(1:Tc));
Template = mean(FirstBeats,2);
%figure, plot(Template);
%pause;

% remoçao de batidas anomalas
index = false(1,m);
for i = Tc+1:m
    Difference = Beats{i} - Template;
    %if norm(Difference) > 1.5
    %    figure, plot(Beats{i});
    %    norm(Difference)
    %    pause;
    %end
    index(i) = norm(Difference) <= 1.5;
end
Beats = cell2mat(Beats(index));
Rpeaks = Rpeaks(index);


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
