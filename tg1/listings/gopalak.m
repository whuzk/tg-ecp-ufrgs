%% Pré-processamento
Signal = ecg_filter_noise_and_baseline(Signal,Fs);
R = ecg_tompkins_rpeak_detection_adapted(Signal,Fs);

%% Extraçao de caracteristicas
C = ecg_extract_hermite_coeff(Signal1,Fs,R,Hncoeff,Hbparam);

%% Classificaçao
IschemicBeats = ecg_classify_ischemic(C);

%% Composiçao das estatisticas
index = ecg_find_close_beats(Assessment.Rpeaks,R);
Diagnosis = Assessment.Diagnosis(index);
Stats = compute_statistics(Diagnosis,IschemicBeats);