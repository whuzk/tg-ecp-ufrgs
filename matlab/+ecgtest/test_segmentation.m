function Result = test_segmentation(Segment, varargin)
import ecgutilities.*

if nargin == 1
    LeadName = 'all';
else
    LeadName = varargin{1};
end

%% count the results
Records = fieldnames(Segment);
count = 0;
for i = 1:numel(Records)
    Record = Segment.(Records{i});
    if strcmp(LeadName,'all')
        Leads = fieldnames(Record);
        for j = 1:numel(Leads)
            count = count + length(Record.(Leads{j}).A);
        end
    elseif isfield(Record,LeadName)
        count = count + length(Record.(LeadName).A);
    end
end

%% group the results
A = false(count,1);
B = false(count,1);
iStart = 1;
iEnd = 0;
for i = 1:numel(Records)
    Record = Segment.(Records{i});
    if strcmp(LeadName,'all')
        Leads = fieldnames(Record);
        for j = 1:numel(Leads)
            Data = Record.(Leads{j});
            iEnd = iEnd + length(Data.A);
            A(iStart:iEnd) = Data.A;
            B(iStart:iEnd) = Data.B;
            iStart = iEnd + 1;
        end
    elseif isfield(Record,LeadName)
        Data = Record.(LeadName);
        iEnd = iEnd + length(Data.A);
        A(iStart:iEnd) = Data.A;
        B(iStart:iEnd) = Data.B;
        iStart = iEnd + 1;
    end
end

%% compute statistics
Result = ecgmath.compute_statistics(A,B);
