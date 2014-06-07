function [Templates,idx] = template_and_artifacts(Beats,Tc)

[M,N] = size(Beats);
Templates = zeros(M,N);
Template = zeros(M,1);
idx = false(1,N);
thr = [0 0];
r = 1/Tc;
figure;
for i = 1:N
    metric = rms(Template - Beats(:,i))
    if (i == 1)
        Template = Beats(:,i);
        plot(Template);
    elseif (i <= Tc || metric < thr(1))
        Template = Template + r * (Beats(:,i) - Template);
        thr = thr + r * ([3*metric 4*metric] - thr);
        plot(Template);
        'very good'
    elseif metric > thr(2)
        idx(i) = true;
        plot(Beats(:,i));
        'bad'
    else
        plot(Beats(:,i));
        'good'
    end
    title(['Beat #' num2str(i)]);
    pause;
    Templates(:,i) = Template;
end