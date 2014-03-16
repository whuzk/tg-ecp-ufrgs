function [Beats,Rpeaks,RR,Template] = preprocess(Signal, Fs)

% detecção dos complexos QRS
tic;
Rpeaks = ecgfilter.detect_qrs(Signal,Fs);
toc;
ecgutilities.plot_signal_r(Signal, Rpeaks);
%pause;

% remoçao de ruido e da linha de base
tic;
[FiltSignal,Delay] = filter_signal(Signal, Rpeaks, Fs);
toc;

% ajuste da localizaçao das batidas
[Radj,Rpeaks,RR] = adjust_peaks(FiltSignal,Rpeaks,Delay);
%ecgutilities.plot_signal_r(FiltSignal, Radj);
%pause;

% extraçao das batidas
Beats = extract_beats(FiltSignal,Radj,RR,2*fix(Fs*0.6)+1);

% construçao de template
Tc = 30;
Template = mean(Beats(:,1:Tc),2);
Beats = Beats(:,Tc+1:end);
Rpeaks = Rpeaks(Tc+1:end);
%ecgutilities.plot_signal(Template, 'Template');
%pause;

% remoçao de batidas anomalas
%
index = detect_artifact_beats(Beats,Template);
Beats = Beats(:,index);
Rpeaks = Rpeaks(index);
%ecgutilities.plot_signal_r(Signal, Rpeaks);
%pause;
%

function [Filtered,Delay] = filter_signal(Signal, Rpeaks, Fs)
[Filtered,T1] = ecgfilter.suppress_noise(Signal, Fs);
%[Filtered,T2] = ecgfilter.suppress_baseline(Filtered, Fs);
[Filtered,T2] = suppress_baseline(Filtered, Rpeaks+T1, Fs);
Delay = T1 + T2;

function [Adjusted,Rpeaks,RR] = adjust_peaks(Signal, Rpeaks, Delay)
RR = diff(Rpeaks(1:end-1));
Rpeaks = Rpeaks(2:end-1);
Adjusted = ecgmath.neighbour_max(Signal,Rpeaks+Delay,5);

function Result = extract_beats(Signal, Rpeaks, RR, FrameSize)
B = ones(1,5)/5;
L = min(RR,filter(B,1,RR));
%figure, plot(Rpeaks, L), grid on;
%figure, plot(Rpeaks, (RR-L).^2), grid on;
%pause
half = fix(L/2);
m = length(Rpeaks);
Result = zeros(FrameSize,m);
for i = 1:m
    left = Rpeaks(i)-half(i);
    right = Rpeaks(i)+half(i);
    Beat = Signal(left:right);
    Result(:,i) = frame_beat(Beat,FrameSize);
end

function Result = detect_artifact_beats(Beats,Template)
m = size(Beats,2);
Result = false(1,m);
for i = 1:m
    Norm = norm(Beats(:,i) - Template);
    Result(i) = Norm <= 4;
    %{
    if Norm > 4
        figure(10), plot(Beats(:,i));
        Norm
        pause;
    end
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
%{
function [Result,Delay] = suppress_noise(Signal, Fs)
Wn = 40 * 2/Fs;
[B,A] = butter(4, Wn);
Result = filter(B, A, Signal);
Delay = fix(mean(grpdelay(B,A)));
%}
function [Result,Delay] = suppress_baseline(Signal, Rpeaks, Fs)
N = length(Signal);
X = Rpeaks - fix(Fs*0.06);
X = unique(X(X > 0));
Y = [0; Signal(X); 0];
Result = Signal - spline(X,Y,(1:N)');
%Result = Signal;
Delay = 0;
