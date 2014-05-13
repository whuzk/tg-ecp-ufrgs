function Result = create_dataset(Beatset, LeadName)
% Criaçao do conjunto de caracteristicas
%
import ecgutilities.*;

%
Records = fieldnames(Beatset);
BeatCount = 0;
for i = 1:numel(Records)
    Record = Beatset.(Records{i});
    if isfield(Record,LeadName)
        Data = Record.(LeadName);
        BeatCount = BeatCount + size(Data.Beats,2);
    end
end
if BeatCount == 0
    error('no beats to be processed');
end

%
Rocha.Features = zeros(BeatCount,14);
Rocha.Diagnosis = zeros(BeatCount,2);
Mohebbi.Features = zeros(BeatCount,20);
Mohebbi.Diagnosis = zeros(BeatCount,2);
Gopalak.Features = zeros(BeatCount,50);
Gopalak.Diagnosis = zeros(BeatCount,2);

%
iStart = 1;
iEnd = 0;
for i = 1:numel(Records)
    Data = Beatset.(Records{i}).(LeadName);
    iEnd = iEnd + size(Data.Beats,2);
    
    % seleciona itens especificos de anotaçao
    STelev = Data.STseg.elevation;
    STdep = Data.STseg.depression;
    Telev = Data.Twave.elevation;
    Tdep = Data.Twave.depression;

    % extraçao de caracteristicas do metodo de Rocha
    Rocha.Features(iStart:iEnd,:) = extract_features('Rocha', Data);
    Rocha.Diagnosis(iStart:iEnd,:) = [STelev, STdep];

    % extraçao de caracteristicas do metodo de Mohebbi
    Mohebbi.Features(iStart:iEnd,:) = extract_features('Mohebbi', Data);
    Mohebbi.Diagnosis(iStart:iEnd,1) = STelev | Telev | STdep | Tdep;
    Mohebbi.Diagnosis(iStart:iEnd,2) = ~(STelev | Telev | STdep | Tdep);

    % extraçao de caracteristicas do metodo de Gopalakrishnan
    Gopalak.Features(iStart:iEnd,:) = extract_features('Gopalak', Data);
    Gopalak.Diagnosis(iStart:iEnd,:) = [STelev | STdep, Telev | Tdep];
    
    iStart = iEnd + 1;
end

%
Result.Rocha = Rocha;
Result.Mohebbi = Mohebbi;
Result.Gopalak = Gopalak;
