function test_template(Beatset)
% Teste dos templates
%
Records = fieldnames(Beatset);
k = numel(Records);
GlobalMeasure = zeros(k,1);
for i = 1:k
    Data = Beatset.(Records{i});
    
    STelev = Data.STseg.elevation;
    STdep = Data.STseg.depression;
    Telev = Data.Twave.elevation;
    Tdep = Data.Twave.depression;
    
    isquemic = STelev | Telev | STdep | Tdep;
    Beats = Data.Beats(:,isquemic);
    m = size(Beats,2);
    if m > 0
        %figure, plot(Data.Template);
        Measure = zeros(m,1);
        for j = 1:m
            NORM = norm(Beats(:,j) - Data.Template);
            if NORM <= 4
                Measure(j) = NORM;
            end
        end
        GlobalMeasure(i) = max(Measure);
    end
end
GlobalMeasure