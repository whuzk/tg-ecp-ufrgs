function [Template,ArtIdx] = detect_artifacts(Beats,Tc,Threshold)

M = size(Beats,2);
if Tc > M
    error('template count be less or equal to number of beats');
end

Template = mean(Beats(:,1:Tc),2);
ArtIdx = false(1,M);
%figure;
for i = Tc+1:M
    ArtIdx(i) = rms(Beats(:,i) - Template) > Threshold;
    %if ArtIdx(i)
    %    plot(Beats(:,i));
    %    rms(Beats(:,i) - Template)
    %    pause;
    %end
end