function Result = jay_point(Beat, Fs)

R = fix(length(Beat)/2);
left = R+fix(0.02*Fs);
right = R+fix(0.12*Fs);
X = Beat(left:right);

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
