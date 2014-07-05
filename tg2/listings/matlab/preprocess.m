function [R,RR,FP,Beats,Template,delay] = preprocess(signal, Fm)
% obtem as amostras e a taxa de amostragem
data = signal.data - signal.inival;
Fs = signal.fs;
% faz a filtragem do sinal
[sigI,delay] = qrs_filter(data, Fs);
[sigD,sigM] = fp_filter(data, Fs, delay);
sigN = noise_filter(data, Fs, Fm, delay);
% faz a deteccao de pontos e de batidas
[R,RR] = detect_qrs(sigI, Fs);
FP = detect_fp(sigD, sigM, R, Fs);
[Beats,FP] = extract_beats(sigN, FP, Fs);
[Templates,idx] = build_template(Beats, 30);
Template = Templates(:,30);
% remove os artefatos
R(idx) = [];
RR(idx) = [];
FP(:,idx) = [];
Beats(:,idx) = [];