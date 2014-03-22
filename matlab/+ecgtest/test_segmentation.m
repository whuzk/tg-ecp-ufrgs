function Result = test_segmentation(Segment)
import ecgutilities.*

%% count the results
Records = fieldnames(Segment);
count = 0;
for i = 1:numel(Records)
    record = Segment.(Records{i});
    count = count + length(record.A);
end

%% group the results
A = false(count,1);
B = false(count,1);
iStart = 1;
iEnd = 0;
for i = 1:numel(Records)
    record = Segment.(Records{i});
    iEnd = iEnd + length(record.A);
    A(iStart:iEnd) = record.A;
    B(iStart:iEnd) = record.B;
    iStart = iEnd + 1;
    ecgmath.compute_statistics(record.A,record.B)
end

%% compute statistics
Result = ecgmath.compute_statistics(A,B);
