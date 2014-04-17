function Result = detect_qrs2(Signal, Fs)

[Filtered,offset] = ecgfilter.sogari_filter(Signal,Fs);
Result = ecgfeatures.sogari_qrs(Filtered,Fs) + offset;
%{
[Result,R2,THs,THn,RRm] = ecgfeatures.sogari_qrs(Filtered,Fs);
Result = Result + offset;
R2 = R2 + offset;
figure;
hold on, grid on;
plot([Filtered THs THn]);
plot(Result, Filtered(Result), 'kx');
plot(R2, Filtered(R2), 'ko');
%figure, plot(RRm);
%}