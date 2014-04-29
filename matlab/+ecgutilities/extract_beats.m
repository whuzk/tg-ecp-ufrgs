function [Beats,R,RR] = extract_beats(data,lead,Fs,FrameSize)
import ecgutilities.*
import ecgfastcode.*
import ecgfilter.*

[R,RR,delay1] = c_prod_detect_qrs_double(data,Fs);
R = adjust_qrs(data,lead,R-floor(delay1),Fs);

[mmd,lap,delay2] = sogari_filter2(data,Fs,50);
R = R + floor(delay2);
F = fiducial_marks(mmd,lap,lead,R,Fs);

data = [zeros(floor(delay2),1); data(1:end-floor(delay2))];
plot_fiducial_marks(data,F);
plot_fiducial_marks(mmd,F);

M = length(F);
Beats = zeros(FrameSize,M);
for i = 1:M
    Beat = data(F(i).P(1):F(i).T(3));
    %figure, plot(Beat);
    Beat = suppress_baseline(Beat,5);
    %figure, plot(Beat);
    r = F(i).R(2) - F(i).P(1) + 1;
    Beat = frame_beat(Beat,r,FrameSize);
    %figure, plot(Beat);
    Beats(:,i) = Beat;
    %pause;
end
R = R - floor(delay2);