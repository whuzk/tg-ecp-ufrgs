function test_segmentation(ECG)
import ecgfilter.*;

Fs = sscanf(ECG.SamplingFrequency, '%d');
Lead = ECG.Signals{1};
Gain = sscanf(Lead.Gain, '%d');
Offset = sscanf(Lead.InitialValue, '%d');
Signal = (Lead.Data - Offset) / Gain;

Signal = suppress_noise(Signal,Fs);
Rpeaks = detect_qrs(Signal,Fs);
plot_signal_r(Signal, Rpeaks);
%length(Rpeaks)
