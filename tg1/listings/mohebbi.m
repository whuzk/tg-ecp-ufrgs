%% Pré-processamento (parte 1)
R = ecg_tompkins_rpeak_detection_original(Signal1,Fs);
R(ecg_extract_artifact_beats(Signal1,Fs,R)) = [];
T = ecg_construct_normal_template(Signal1,R,Roffset,Tsize);

%% Pré-processamento (parte 2)
R = ecg_tompkins_rpeak_detection_original(Signal2,Fs);
R(ecg_extract_artifact_beats_based_on_template(Signal2,R,Roffset,T,Tthresh)) = [];
Signal2 = ecg_baseline_wander_removal(Signal2,Fs,R);
J = ecg_detect_jay_points(Signal2,Fs,R,Javgw);

%% Preparaçao dos dados
RT = Tsize/2-Roffset+1;
JT = ecg_detect_jay_points(T,Fs,RT,Javgw);
TemplateST = T(JT:min(length(T),JT+STlen-1));
STdata = ecg_prepare_st_segments(Signal2,J,TemplateST);

%% Classificaçao
IschemicBeats = ecg_classify_ischemic_beats(STdata);

%% Composiçao das estatisticas
index = ecg_find_close_beats(Assessment.Rpeaks,R);
Diagnosis = Assessment.Diagnosis(index);
Stats = compute_statistics(Diagnosis,IschemicBeats);