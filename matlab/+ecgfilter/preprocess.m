function [Beats,Rpeaks,RR,Template] = preprocess(Signal, Fs)

% detecção dos complexos QRS
Rpeaks = ecgfeatures.detect_qrs(Signal,Fs);

% remoçao de ruido
Signal = ecgfilter.suppress_noise(Signal,Fs);

% ajuste dos picos
Rpeaks = ecgmath.neighbour_max(Signal,Rpeaks,5);
RR = diff([0; Rpeaks]);

% extraçao das batidas
FrameSize = 2*fix(Fs*0.6)+1;
Beats = ecgutilities.extract_beats(Signal,Rpeaks,RR,FrameSize);
%ecgutilities.plot_signal(Beats(:,1), 'Beat #1');
%pause;

% remoçao da linha de base
Beats = remove_baseline(Beats);
%ecgutilities.plot_signal(Beats(:,1), 'Beat #1');
%pause;

% construçao do template
TemplateCount = 30;
Template = mean(Beats(:,1:TemplateCount),2);
%ecgutilities.plot_signal(Template, 'Template');
%pause;

% remoçao de batidas anomalas
Threshold = 4;
index = detect_artifact_beats(Beats,Template,Threshold);
Beats = Beats(:,index);
Rpeaks = Rpeaks(index);
%ecgutilities.plot_signal_r(Signal, Rpeaks);
%pause;


function Beats = remove_baseline(Beats)
for i = 1:size(Beats,2)
    Beats(:,i) = ecgfilter.suppress_baseline(Beats(:,i),5);
end

function Result = detect_artifact_beats(Beats,Template,Threshold)
m = size(Beats,2);
Result = false(1,m);
for i = 1:m
    Result(i) = norm(Beats(:,i) - Template) <= Threshold;
end