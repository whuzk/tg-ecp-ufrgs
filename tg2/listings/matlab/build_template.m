function [Templates,ArtifactIndices] = build_template(Beats,Tc)
[M,N] = size(Beats);
Templates = zeros(M,N);
Template = zeros(M,1);
ArtifactIndices = false(1,N);
thr = [0 0];
r = 1/Tc;
for i = 1:N
    metric = rms(Template - Beats(:,i))
    if (i == 1)
        % primeiro batimento
        Template = Beats(:,i);
    elseif (i <= Tc || metric < thr(1))
        % batimento excelente
        Template = Template + r * (Beats(:,i) - Template);
        thr = thr + r * ([3*metric 4*metric] - thr);
    elseif metric > thr(2)
        % batimento ruim
        ArtifactIndices(i) = true;
    else
        % batimento bom
    end
    Templates(:,i) = Template;
end