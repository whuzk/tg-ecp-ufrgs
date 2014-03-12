function [Beats,Rpeaks,RR,Template] = ecg_preprocess(Signal, Fs)

% detecção dos complexos QRS
tic;
Rpeaks = ecgutilities.ecg_detect_qrs(Signal,Fs);
toc;
%ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;

% remoçao de ruido e da linha de base
tic;
[Signal,T1] = suppress_noise(Signal,Fs);
[Signal,T2] = suppress_baseline(Signal,Rpeaks,Fs);
toc;

% ajuste da localizaçao das batidas
[Rpeaks,RR] = adjust_peaks(Signal,Rpeaks,T1+T2);
%ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;

% extraçao das batidas
Beats = extract_beats(Signal,Rpeaks,RR,2*fix(Fs*0.6)+1);

% construçao de template
Tc = 30;
Template = mean(Beats(:,1:Tc),2);
Beats = Beats(:,Tc+1:end);
Rpeaks = Rpeaks(Tc+1:end);
%ecgutilities.ecg_plot(Template, 'Template');
%ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;

% remoçao de batidas anomalas
%{
tic;
index = detect_artifact_beats(Beats,Template);
toc;
Beats = Beats(:,index);
Rpeaks = Rpeaks(index);
ecgutilities.ecg_plot_r(Signal, Rpeaks);
%pause;
%}

function [Result,Delay] = suppress_noise(Signal, Fs)
Wn = 40 * 2/Fs;
[B,A] = butter(4, Wn);
Result = filter(B, A, Signal);
Delay = fix(mean(grpdelay(B,A)));

function [Result,Delay] = suppress_baseline(Signal, Rpeak, Fs)
N = length(Signal);
X = Rpeak - fix(Fs*0.06);
Y = [0; Signal(X); 0];
Result = Signal - spline(X,Y,(1:N)');
%Result = Signal;
Delay = 0;

function [Rpeaks,RR] = adjust_peaks(Signal,Rpeaks,Delay)
N = length(Signal);
Rpeaks = Rpeaks + Delay;
Rpeaks(Rpeaks > N) = [];
Rpeaks = ecgmath.neighbour_max(Signal,Rpeaks,5);
RR = diff(Rpeaks(1:end-1));
Rpeaks = Rpeaks(2:end-1);

function Result = extract_beats(Signal, Rpeaks, RR, MaxLength)
half = fix(min(RR,MaxLength)/2);
m = length(Rpeaks);
Result = zeros(MaxLength,m);
for i = 1:m
    left = Rpeaks(i)-half(i);
    right = Rpeaks(i)+half(i);
    Beat = Signal(left:right);
    Result(:,i) = frame_beat(Beat,MaxLength);
end

function Result = detect_artifact_beats(Beats,Template)
m = size(Beats,2);
Result = false(1,m);
for i = 1:m
    Norm = norm(Beats(:,i) - Template);
    %{
    if Norm > 3
        figure(10), plot(Beats(:,i));
        Norm
        pause;
    end
    %}
    Result(i) = Norm <= 3;
    %{
    Correlation = corr(Beats(:,i), Template);
    if Correlation < 0.75
        figure(10), plot(Beats(:,i));
        Correlation
        pause;
    end
    Result(i) = Correlation >= 0.75;
    %}
end

function Result = frame_beat(Signal, l)
D = l-length(Signal);
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    Result = Signal(1+a:end-b);
else
    Result = [zeros(a,1); Signal; zeros(b,1)];
end
