function stats = generate_sel_stats(basedir)
stats = cell(5,1);
for i = 1:5
    load([basedir 'superset' num2str(i) '.mat']);
    stats{i}.RochaST = create_stat_table(RochaSet);
    stats{i}.MohebbiST = create_stat_table(MohebbiSet);
    stats{i}.GopalakST = create_stat_table(GopalakSTSet);
    stats{i}.GopalakT = create_stat_table(GopalakTSet);
    clear RochaSet MohebbiSet GopalakSTSet GopalakTSet
end

function Result = create_stat_table(datasets)
fields = fieldnames(datasets);
stats = zeros(length(fields),3);
for i = 1:length(fields)
    t = datasets.(fields{i}).targets;
    stats(i,1) = length(t);
    stats(i,2) = length(find(t < 0));
    stats(i,3) = length(find(t > 0));
end
Result = [...
    table(fields, 'VariableNames', {'LeadName'})...
    array2table(stats, 'VariableNames', ...
    {'TotalCount' 'NormalCount' 'IschemicCount'})
];
ratio = Result.NormalCount ./ Result.TotalCount;
Result.NormalPercent = num2str(ratio*100,'%.1f');
ratio = Result.IschemicCount ./ Result.TotalCount;
Result.IschemicPercent = num2str(ratio*100,'%.1f');