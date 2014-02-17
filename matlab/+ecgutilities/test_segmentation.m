function test_segmentation(ECG)
import ecgutilities.*;

Fs = sscanf(ECG.SamplingFrequency, '%d');
Lead = ECG.Signals{1};
Gain = sscanf(Lead.Gain, '%d');
Offset = sscanf(Lead.InitialValue, '%d');
Signal = (Lead.Data - Offset) / Gain;

Signal = suppress_noise(Signal,Fs);
Rpeaks = ecg_segment(Signal,Fs);
ecg_plot_r(Signal, Rpeaks);
%length(Rpeaks)


function Result = suppress_noise(Signal, Fs)
Wn = 40 * 2/Fs;
[B,A] = butter(4, Wn);
Result = filter(B, A, Signal);
