function [Beats,Rd,RR,Template] = preprocess(Signal)

% inicializacao
data = Signal.data - Signal.inival;

% detecção dos complexos QRS
[~,Rd,~,~,~,~,delay1] = ecgfastcode.detect_qrs_double(data,Signal.fs);

% remoçao de ruido
[data,delay2] = ecgfilter.suppress_noise(data,Signal.fs);

% atraso do sinal e ajuste dos picos
[data,Rd] = adjust(data,Rd,floor(delay1-delay2),Signal.lead);
%ecgutilities.plot_signal_r(data, Rd);
%pause;

% extraçao das batidas e remoçao da linha de base
RR = diff([0; Rd]);
FrameSize = 2*fix(Signal.fs*0.6)+1;
Beats = ecgutilities.extract_beats(data,Rd,RR,FrameSize);
%ecgutilities.plot_signal(Beats(:,1), 'Beat #1');
%pause;

% construçao do template
TemplateCount = 30;
Template = mean(Beats(:,1:TemplateCount),2);
%ecgutilities.plot_signal(Template, 'Template');
%pause;

% remoçao de batidas anomalas
index = detect_artifact_beats(Beats,Template,20);
%ecgutilities.plot_signal_r(data, Rd(~index));
%pause;

% salva o resultado
Beats = Beats(:,~index);
Rd = Rd(~index)-floor(delay1);


function Result = detect_artifact_beats(Beats,Template,Threshold)
m = size(Beats,2);
Result = false(1,m);
%figure;
for i = 1:m
    Result(i) = rms(Beats(:,i) - Template) > Threshold;
    %if Result(i)
    %    plot(Beats(:,i));
    %    rms(Beats(:,i) - Template)
    %    pause;
    %end
end

function [data,Rd] = adjust(data,Rd,delay,lead)
data = [zeros(delay,1); data(1:end-delay)];
if ismember(lead,{'MLIII'})
    Rd = ecgmath.neighbour_max(-data,Rd,5);
else
    Rd = ecgmath.neighbour_max(data,Rd,5);
end