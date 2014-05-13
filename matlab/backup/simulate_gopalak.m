function simulate_gopalak

import utilities.*;
global EDB Hs;
close all;

%% obtem informaçoes do ECG
Records = fieldnames(EDB);
ECG = EDB.(Records{4});
[Fs,Bp,~,~,Data] = ecg_interpret(ECG);
Hs = 2*fix(Fs*0.6)+1;

j = 2;
SignalID = num2str(j-1);
Signal = Data(1:10000,j);
Rpeaks = Bp(Bp < 10000);

% obtem os diagnosticos
ST = ecg_get_diagnosis(ECG.Annotations, Rpeaks, 's', 'ST', SignalID);
T = ecg_get_diagnosis(ECG.Annotations, Rpeaks, 'T', 'T', SignalID);
Diagnosis = ST.elevation | ST.depression | T.elevation | T.depression;

%% filtros
Wn = 40 * 2/Fs;
%Wn = [0.5 40] * 2/Fs;
[B1,A1] = butter(4, Wn);

% low-pass (5T delay)
A = [1 -2 1];
B([1 7 13]) = 1/36*[1 -2 1];
T1 = tf(B,A);
clear A B;

% high-pass (16T delay)
A = [1 -1];
B([1 17 18 33]) = [-1/32 1 -1 1/32];
T2 = tf(B,A);
clear A B;

% derivative (2T delay)
A = 1;
B = 0.1*[2 1 0 -1 -2];
T3 = tf(B,A);
clear A B;

% cascade (23T delay)
T = T1*T2*T3;
A2 = fliplr(T.den{1});
B2 = T.num{1};

% integration (18T delay)
Ws = fix(0.15*Fs);
B = ones(1,Ws)/Ws;
T4 = tf(B,1);

% smoothing (1T delay)
B = ones(1,3)/3;
T5 = tf(B,1);

% cascade (19T delay)
T = T4*T5;
A3 = fliplr(T.den{1});
B3 = T.num{1};

% hermite
H = ecgmath.hermite(Hs, 50, 1.5);

%% objetos da simulaçao
N = 1000;
Signal11 = zeros(N,1);
Signal12 = zeros(N,1);
Signal13 = zeros(N,1);
Aux1 = zeros(N,1);
Signal21 = zeros(Hs,1);
Signal22 = zeros(Hs,1);
Signal23 = zeros(Hs,1);
Aux2 = zeros(Hs,1);

figure(1);

% grafico 1
subplot(3,1,1); grid on;
title('Original signal');
lh11 = line((1:N)',Signal11);
lhD = line(0,0, 'marker','o', 'markersize',40, 'linestyle','none', 'color','black');

% grafico 2
subplot(3,1,2); grid on;
title('Filtered signal');
lh12 = line((1:N)',Signal12);
lhR = line(0,0, 'marker','.', 'markersize', 5, 'linestyle','none', 'color','black');
lhI = line(0,0, 'marker','o', 'markersize',10, 'linestyle','none', 'color','black');

% grafico 3
subplot(3,1,3); grid on;
title('Segmentation');
lh13 = line((1:N)',Signal13);
lhIT = line([1; N],[0; 0], 'linestyle','--', 'color','green');
lhIR = line(0,0, 'marker','x', 'markersize',10, 'linestyle','none', 'color','red');

figure(2);

% grafico 1
subplot(3,1,1); grid on;
title('Cardiac beat');
lh21 = line(0,0);
lhB = line(0,0, 'LineWidth',2, 'color','red');

% grafico 2
subplot(3,1,2); grid on;
title('Gradient');
lh22 = line(0,0);

% grafico 3
subplot(3,1,3); grid on;
title('Hermite expansion');
xlim([1 Hs]);
lh23 = line((1:Hs)',Signal23);
lhH = line(1:Hs,Signal23, 'marker','.', 'linestyle','none', 'color','black');


%
%% inicio
step = 1/Fs;
avg1 = 0;
avg2 = 0;
NumBeats = 0;
Stresh = 0;
Ntresh = 0;
Tresh1 = 0;
skip = false;
Fdelay = 23+19;
R = zeros(0,1);
RR = zeros(0,1);
R2 = zeros(0,1);
I = false(0,1);
D = zeros(0,1);
F = zeros(10,50);

for i = 1:length(Signal)
    % um passo do processamento
    s1 = tic;
    Signal11 = [Signal11(2:end); Signal(i)];
    Signal12 = [Signal12(2:end); apply_filter(Signal11, Signal12, B1, A1)];
    Aux1 = [Aux1(2:end); apply_filter(Signal11, Aux1, B2, A2)];
    Signal13 = [Signal13(2:end); apply_filter(Aux1.^2, Signal13, B3, A3)];
    [Stresh,Ntresh,Tresh1,skip,ok] = ...
        update_detection(Signal13, Stresh, Ntresh, Tresh1, skip);
    [R,R2,I] = update_Rpoints(Signal12, R, R2, I, Fdelay, ok);
    if ok
        RR = diff(R);
    end
    
    t2 = 0;
    if ~isempty(RR)
        L = min(Hs,RR(end));
        half = fix(L/2);
        if R(end) == N-half
            Signal21 = Signal12(end-2*half:end);
            %Signal21 = process_beat(Signal21, B1, A1);
            F = [F(2:end,:); zeros(1,50)];
            s2 = tic;
            [Signal22,Aux2,Signal23,F(end,:),I(end)] = ...
                OnNewBeat(Signal21, Fs, half+1, H, F);
            t2 = toc(s2);
            avg2 = avg2 + t2/L;
            NumBeats = NumBeats + 1;
            pause;
        end
    end
    avg1 = avg1 + toc(s1) - t2;
    
    % diagnosticos (nao conta no tempo)
    D = update_diagnosis(Signal11, Rpeaks, Diagnosis, i, D);
    
    % atualiza o grafico (nao conta no tempo)
    %if mod(i,2) == 0
        set(lh11, 'ydata', Signal11);
        set(lhD,  'xdata', D, 'ydata', Signal11(D));
        set(lh12, 'ydata', Signal12);
        set(lhR, 'xdata', R, 'ydata', Signal12(R));
        set(lhI, 'xdata', R(I), 'ydata', Signal12(R(I)));
        set(lh13, 'ydata', Signal13);
        set(lhIT, 'ydata', [Tresh1,Tresh1]);
        set(lhIR, 'xdata', R2, 'ydata', Signal13(R2));
        set(lh21, 'xdata', 1:length(Signal21), 'ydata', Signal21);
        set(lhB, 'xdata', 1:length(Aux2), 'ydata', Aux2);
        set(lh22, 'xdata', 1:length(Signal22), 'ydata', Signal22);
        set(lh23, 'ydata', Signal23);
        set(lhH, 'ydata', H*F(end,:)');
        drawnow;
    %end
end

avg1 = avg1 / N;
avg2 = avg2 / NumBeats;
avg = avg1 + avg2;
disp(['sampling period of ECG: ' num2str(step*1000) ' ms']);
disp(['average time per-sample: ' num2str(avg*1000) ' ms']);
disp(['max. freq. supported: ' num2str(fix(1/avg)) ' Hz']);


function [Gradient,Y,X,newF,newI] = OnNewBeat(Beat, Fs, R, H, F)
global Hs;

[i,a,Gradient] = utilities.ecg_baseline_points(Beat, Fs, R);

n = length(Beat);
Y = spline(i,a,(1:n)');
X = Beat - Y;

D = Hs-n;
a = fix(abs(D)/2);
b = abs(D)-a;
if D <= 0
    X = X(1+a:end-b);
else
    X = [zeros(a,1); X; zeros(b,1)];
end

newF = X'*H;
B = ones(1,10)/10;
temp = [F(2:end,:); newF];
newI = classify_beat(B*temp);

%{
function Result = process_beat(Beat, B1, A1)

Temp = filter(B1, A1, Beat(end:-1:1));
Temp = Temp(end:-1:1);
pks = -findpeaks(-Temp);
Result = Temp-mean(pks);
%}

function Result = classify_beat(F)
global GopalakNets;

k = length(GopalakNets);
decision = 0;
for i = 1:k
    O = sim(GopalakNets{i}, F')';
    C = O(:,1) > 0 | O(:,2) > 0;
    decision = decision + double(C);
end
Result = decision > k/2;


function [R,R2,I] = update_Rpoints(Signal, R, R2, I, delay, HasNewR)

% atualiza os picos R
R = R - 1;
R2 = R2 - 1;
if ~isempty(R) && R(1) == 0
    R(1) = [];
    R2(1) = [];
    I(1) = [];
end

% novo pico R
if HasNewR
    pos = length(Signal)-1;
    x = Signal(pos-2*delay:pos);
    [~,j] = max(x);
    R = [R; pos-2*delay+j-1];
    R2 = [R2; pos];
    I = [I; false];
end


function D = update_diagnosis(Signal, Rpeaks, Diagnosis, i, D)

N = length(Signal);
D = D - 1;
if ~isempty(D) && D(1) == 0
    D(1) = [];
end
index = find(i == Rpeaks, 1, 'first');
if ~isempty(index) && Diagnosis(index)
    D = [D; N];
end


function [Stresh,Ntresh,Tresh1,skip,ok] = update_detection(X, Stresh, Ntresh, Tresh1, skip)

i = length(X)-1;
ok = false;
if skip
    if X(i) < Tresh1
        skip = false;
    end
else
    a = X(i-1);
    b = X(i);
    c = X(i+1);
    if (b-a) > 0 && (c-b) <= 0 && abs(c-2*b+a) > 1E-6
        PEAK = X(i);
        if PEAK > Tresh1
            Stresh = 0.25*PEAK + 0.75*Stresh;
            skip = true;
            ok = true;
        else
            Ntresh = 0.25*PEAK + 0.75*Ntresh;
        end
        Tresh1 = Ntresh + 0.25*(Stresh-Ntresh);
    end
end

function Result = apply_filter(X, Y, B, A)

na = length(A);
nb = length(B);
Result = (B*X(end:-1:end-nb+1) - A(2:end)*Y(end:-1:end-na+2))/A(1);
