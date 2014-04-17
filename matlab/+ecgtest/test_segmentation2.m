function Result = test_segmentation2(Segment, LeadName)
import ecgutilities.*

%% count the results
Records = fieldnames(Segment);
count = 0;
for i = 1:numel(Records)
    Record = Segment.(Records{i});
    Leads = fieldnames(Record);
    for j = 1:numel(Leads)
        if ~strcmp(Leads{j},LeadName)
            count = count + length(Record.(Leads{j}).A);
        end
    end
end

%% group the results
A = false(count,1);
B = false(count,1);
iStart = 1;
iEnd = 0;
for i = 1:numel(Records)
    Record = Segment.(Records{i});
    Leads = fieldnames(Record);
    for j = 1:numel(Leads)
        if ~strcmp(Leads{j},LeadName)
            Data = Record.(Leads{j});
            iEnd = iEnd + length(Data.A);
            A(iStart:iEnd) = Data.A;
            B(iStart:iEnd) = Data.B;
            iStart = iEnd + 1;
        end
    end
end

%% compute statistics
Result = ecgmath.compute_statistics(A,B);
