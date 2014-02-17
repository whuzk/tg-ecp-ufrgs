function simulate

% Simluaçao de tempo-real dos metodos de detecçao
import utilities.*;
global EDB;
close all;

% obtem informaçoes do ECG
names = fieldnames(EDB);
RecordName = names{4};
ECG = EDB.(RecordName);
Fs = sscanf(ECG.SamplingFrequency, '%d');
BP = str2double(ECG.Annotations.Sample);

% obtem os dados da derivaçao
j = 2;
Lead = ECG.Signals{j};
Gain = sscanf(Lead.Gain, '%d');
Offset = sscanf(Lead.InitialValue, '%d');
Data = (Lead.Data - Offset) / Gain;
SignalID = num2str(j-1);

%
%Data = Data(73001:end);
%Rpeaks = BP(73000 < BP);
%Data = Data(70001:end);
%Rpeaks = BP(70000 < BP);
%Data = Data(1:3000);
%Rpeaks = BP(BP < 3000);
%Data = Data(200001:end);
%Rpeaks = BP(200000 < BP);
Data = Data(223001:end);
Rpeaks = BP(223000 < BP);
%Rpeaks = BP;
       
% obtem os diagnosticos
ST = ecg_get_diagnosis(ECG.Annotations, Rpeaks, 's', 'ST', SignalID);
T = ecg_get_diagnosis(ECG.Annotations, Rpeaks, 'T', 'T', SignalID);
Diagnosis = ST.elevation | ST.depression | T.elevation | T.depression;

%Rpeaks = Rpeaks - 73000;
%Rpeaks = Rpeaks - 70000;
%Rpeaks = Rpeaks - 200000;
Rpeaks = Rpeaks - 223000;

%% Simulaçao
N = 1000;
Signal1 = zeros(N,1);
Signal2 = zeros(N,1);

% grafico 1
subplot(2,1,1); grid on;
title('Original signal');
lh1 = line((1:N)',Signal1);
lhD = line(0,0, 'marker','o', 'markersize',40, 'linestyle','none', 'color','black');

% grafico 2
subplot(2,1,2); grid on;
title('Output signal');
lh2 = line((1:N)',Signal2);
lhR = line(0,0, 'marker','.', 'markersize',5, 'linestyle','none', 'color','black');
lhI = line(0,0, 'marker','o', 'markersize',10, 'linestyle','none', 'color','black');
lhB = line([0; 0],[-1; 2], 'linestyle','--', 'color','red');

% filtros
Wn = [0.5 40] * 2/Fs;
[B,A] = butter(4, Wn);
B2 = ones(1,10)/10;
B3 = ones(1,50)/50;

% inicio
%nf = 14;
%nf = 20;
nf = 50;
step = 1/Fs;
delay = 0;
begin = N;
R = zeros(0,1);
I = false(0,1);
D = zeros(0,1);
M = zeros(50,1);
F = zeros(10,nf);
%{
ok = false;
count = 0;
Roff = fix(0.12*Fs);
STlen = fix(0.16*Fs);
Template = zeros(201,1);
%}
BL = 0;
for i = 1:length(Data)
    % um passo do processamento
    tic;
    [Signal1,Signal2,Out] = update_signals(Signal1, Signal2, Data(i), B, A);
    [begin,R,I,D,BL,new] = update_points(Out, Fs, Rpeaks, Diagnosis, i, begin, R, I, D, M, BL, B3);
    Out = Out - BL;
    %{
    if i > 10*Fs && ~ok
        ok = true;
        Template = Template/count;
        Rtemplate = 100-Roff+1;
        Jtemplate = ecgmohebbi.ecg_detect_jay_points(Template, Fs, Rtemplate);
        TemplateST = Template(Jtemplate:min(201,Jtemplate+STlen-1));
        STmodel = (TemplateST(1:2:end) + TemplateST(2:2:end))/2;
        disp('Template ok');
        %{
        figure;
        hold on; grid on;
        plot(Template);
        plot(Rtemplate, Template(Rtemplate), 'kx');
        plot(Jtemplate, Template(Jtemplate), 'kx');
        title('Template');
        ecg_plot(TemplateST, 'Template ST segment');
        pause;
        %}
    end
    %}
    if new && size(R,1) > 1
        %[F,I] = update_rocha_caracteristics(Out,Fs,R(end-1:end),F,I,B2);
        %{
        if ~ok
            [Template,count] = update_mohebbi_template(Out,R+Roff,Template,count);
        else
            [F,I] = update_mohebbi_caracteristics(Out,Fs,R(end-1:end),F,I,B2,Template,STmodel,Roff,STlen);
        end
        %}
        [F,I] = update_gopalak_caracteristics(Out,Fs,R(end-1:end),F,I,B2);
    end
    delay = delay + max(0,toc-step);
    
    % atualiza o grafico
    subplot(2,1,1);
    set(lh1, 'ydata', Signal1);
    set(lhD, 'xdata', D, 'ydata', Signal1(D));
    subplot(2,1,2);
    set(lh2, 'ydata', Out);
    set(lhR, 'xdata', R, 'ydata', Out(R));
    set(lhI, 'xdata', R(I), 'ydata', Out(R(I)));
    set(lhB, 'xdata', [begin,begin]);
    drawnow;
end

if delay > 0
    warning(['total delay: ' num2str(delay)]);
end


function [X, Y, Result] = update_signals(X, Y, NewValue, B, A)

X = [X(2:end); NewValue];
Y = [Y(2:end); 0];
na = length(A);
nb = length(B);
Y(end) = (B*X(end:-1:end-nb+1) - A(2:end)*Y(end-1:-1:end-na+1))/A(1);

Result = filter(B, A, Y(end:-1:1));
Result = Result(end:-1:1);

%temp = filter(B, A, Y(end:-1:400));
%Result = [Y(1:399); temp(end:-1:1)];


function [begin,R,I,D,BL,new] = update_points(Signal, Fs, Rpeaks, Diagnosis, i, begin, R, I, D, M, BL, B)
import ecggopalak.*;
import utilities.*;

N = length(Signal);

% atualiza os minimos
a = Signal(end-20);
b = Signal(end-19);
c = Signal(end-18);
if (b-a) < 0 && (c-b) >= 0 && abs(c-2*b+a) > 1E-4
    M = [M(2:end); Signal(N-19)];
    BL = B*M(end:-1:1);
end
%
% atualiza os diagnosticos
D = D - 1;
if ~isempty(D) && D(1) == 0
    D(1) = [];
end
index = find(i == Rpeaks, 1, 'first');
if ~isempty(index) && Diagnosis(index)
    D = [D; N];
end
%
% atualiza os picos de onda R
R = R - 1;
if ~isempty(R) && R(1) == 0
    R(1) = [];
    I(1) = [];
end

% dedecta um novo pico de onda R
if begin == N-2*Fs
    newR = ecg_segment(Signal(begin:end), Fs, 1);
    newR = newR + begin - 1;
else
    begin = begin - 1;
    newR = [];
end

% verifica se foi detectada uma batida nova
if ~isempty(newR)
    R = [R; newR(1)];
    I = [I; false];
    begin = min(N,newR(1)+fix(0.1*Fs));
    new = true;
else
    new = false;
end


function [F,I] = update_rocha_caracteristics(Signal,Fs,R,F,I,B)
import ecgrocha.*;

RR = diff(R);
L = min(Fs,RR);
h = fix(L/2);
j = (-h:h)+R(end);
X = Signal(j);

[A,S,J] = ecg_extract_st_deviation(X, Fs, h+1, RR);
[C1,C2] = ecg_hermite_coefficients(X, h+1, S, J, L);
B1 = X(A) - X(S);
B2 = X(J) - X(S);
%ecg_plot_st_deviation(X, A, S, J);
%ecg_plot_hermite_expansion(X, h+1, S, J, L, C1, C2);
%pause;

newF = [B1 B2 C1 C2];
F = [F(2:end,:); newF];
newF = B*F(end:-1:1,:);
I(end) = ecg_classify_ischemic_beats(newF);


function [Template,count] = update_mohebbi_template(Signal,R,Template,count)

L = length(Template);
h = fix(L/2);
j = (-h:h)+R(end);
X = Signal(j);
Template = Template + X;
count = count + 1;


function [F,I] = update_mohebbi_caracteristics(Signal,Fs,R,F,I,B,Template,STmodel,Roff,STlen)
import ecgmohebbi.*;
import utilities.*;

L = length(Template);
h = fix(L/2);
j = (-h:h)+R(end)+Roff;
X = Signal(j);

Difference = X-Template;
if norm(Difference) > 1.0
    return;
end

RR = diff(R);
L = min(Fs,RR);
h = fix(L/2);
j = (-h:h)+R(end);
X = Signal(j);

J = ecg_detect_jay_points(X, Fs, h+1);
if J == 0 || J > length(X)-STlen+1
    return
end

N = length(X);
l = 2*length(STmodel);
ST = X(J:min(N,J+l-1));
ST = (ST(1:2:end) + ST(2:2:end))/2;
%{
figure;
hold on; grid on;
plot(X);
plot(h+1, X(h+1), 'kx');
plot(J, X(J), 'kx');
title('R and J point');
ecg_plot(ST, 'ST Segment');
pause;
%}
newF = (STmodel - ST)';
F = [F(2:end,:); newF];
newF = B*F(end:-1:1,:);
I(end) = ecg_classify_ischemic_beats(newF);


function [F,I] = update_gopalak_caracteristics(Signal,Fs,R,F,I,B)
import ecggopalak.*;
import ecgmath.*;

RR = diff(R);
L = min(Fs,RR);
h = fix(L/2);
j = (-h:h)+R(end);
X = Signal(j);

newF = X'*hermite(length(X), 50, 1.5);
%ecg_plot_hermite_expansion(X, newF, h+1, L);
%pause;

F = [F(2:end,:); newF];
newF = B*F(end:-1:1,:);
I(end) = ecg_classify_ischemic_beats(newF);
