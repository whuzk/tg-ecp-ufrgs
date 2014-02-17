function Result = detect_jay_point(Signal, Fs)

R = fix(length(Signal)/2);
left = R+fix(0.02*Fs);
right = R+fix(0.12*Fs);
X = Signal(left:right);

lf = fix(0.02*Fs);
B = ones(1,lf)/lf;
avg = filtfilt(B,1,X);
grad = gradient(avg);
i = find(grad <= 2.5/Fs, 1, 'first');

if ~isempty(i)
    Result = i+left-1;
else
    Result = left;
end
