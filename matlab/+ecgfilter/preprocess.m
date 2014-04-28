function [Beats,Rd,RR,Template] = preprocess(Signal)

% inicializacao
data = Signal.data - Signal.inival;

% detec��o dos complexos QRS
[~,Rd,~,~,~,~,delay1] = ecgfastcode.detect_qrs_double(data,Signal.fs);

% remo�ao de ruido
[data,delay2] = ecgfilter.suppress_noise(data,Signal.fs);

% atraso do sinal e ajuste dos picos
Rd = adjust_qrs(data,Rd,delay1-delay2,Signal.lead);
%ecgutilities.plot_signal_r(data, Rd);
%pause;

% calculo dos intervalos RR
RR = diff(Rd);

% extra�ao das batidas e remo�ao da linha de base
FrameSize = 2*fix(Signal.fs*0.6)+1;
Beats = ecgutilities.extract_beats(data,Rd,FrameSize);
%ecgutilities.plot_signal(Beats(:,1), 'Beat #1');
%pause;

% constru�ao do template
TemplateCount = 30;
Template = mean(Beats(:,1:TemplateCount),2);
%ecgutilities.plot_signal(Template, 'Template');
%pause;

% remo�ao de batidas anomalas
index = detect_artifact_beats(Beats,Template,20);
%ecgutilities.plot_signal_r(data, Rd(~index));
%pause;

% salva o resultado
RR = RR(~index);
Beats = Beats(:,~index);
Rd = Rd(~index)-floor(delay2);


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

function Rd = adjust_qrs(data,Rd,delay,lead)
Rd = Rd - floor(delay);
if ismember(lead,{'MLIII'})
    Rd = ecgmath.neighbour_max(-data,Rd,5);
else
    Rd = ecgmath.neighbour_max(data,Rd,5);
end