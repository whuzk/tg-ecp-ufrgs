function [Beats,Rd,RR,Template] = preprocess(Signal)
import ecgutilities.*
import ecgfastcode.*
import ecgfilter.*

% inicializacao
data = Signal.data - Signal.inival;

% detecção dos complexos QRS
[Rd,RR,delay1] = prod_detect_qrs_double(data,Signal.fs);

% remoçao de ruido
[data,delay2] = suppress_noise(data,Signal.fs);

% ajuste dos picos
Rd = Rd - floor(delay1 - delay2);
Rd = adjust_qrs(data,Signal.lead,Rd,Signal.fs);

% extraçao das batidas e remoçao da linha de base
FrameSize = 2*floor(Signal.fs*0.6)+1;
Beats = extract_beats(data,Signal.lead,Signal.fs,Rd,RR,FrameSize);

% construçao de template e remoçao de batidas anomalas
[Template,index] = detect_artifacs(Beats,30,30);
%{
ecgutilities.plot_signal_r(data, Rd);
ecgutilities.plot_signal(Beats(:,1), 'Beat #1');
ecgutilities.plot_signal(Template, 'Template');
ecgutilities.plot_signal_r(data, Rd(~index));
%}
% salva o resultado
RR = RR(~index);
Beats = Beats(:,~index);
Rd = Rd(~index)-floor(delay2);
