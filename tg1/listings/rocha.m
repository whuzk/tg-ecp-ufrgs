%% Pré-processamento
Signal = ecg_filter_noise(Signal,Fs);
[Beats,P,Q,R,S,T] = ecg_sun_segment_adapted(Signal,Fs,Scale);
PVC = ecg_extract_pvc_beats(Signal,LeadName,Fs,P,Q,R,S,T);
Beats(PVC,:) = [];
P(PVC,:) = [];
Q(PVC,:) = [];
R(PVC,:) = [];
S(PVC,:) = [];
T(PVC,:) = [];
Signal = ecg_remove_baseline(Signal,Beats);

%% Extração de características
[A,B,J] = ecg_extract_st_deviation(Signal,Fs,R,Beats);
[C1,C2] = ecg_hermite_coefficients(Signal,Fs,P,J,T);
F1 = [A B];
F2(:,:,1) = C1;
F2(:,:,2) = C2;

%% Classificação
I1 = ecg_ischemic_st_elevation(F1,F2);
I2 = ecg_ischemic_st_depression(F1,F2);
IschemicBeats = I1 | I2;
EP = ecg_ischemic_episode_detection(IschemicBeats,Isize);

%% Composiçao das estatisticas
index = ecg_find_close_beats(Assessment.Rpeaks,R(:,2));
Diagnosis = Assessment.Diagnosis(index);
Stats = compute_statistics(Diagnosis,IschemicBeats);