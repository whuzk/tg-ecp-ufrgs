function tables = generate_col_stats(basedir)

tables = cell(20,1);
for i = 1:5
    disp(['Getting statistics for configuration ' num2str(i)]);
    load([basedir 'networks' num2str(i) '.mat']);
    tables{4*(i-1)+1} = get_stat_table(RochaNet);
    tables{4*(i-1)+2} = get_stat_table(MohebbiNet);
    tables{4*(i-1)+3} = get_stat_table(GopalakSTNet);
    tables{4*(i-1)+4} = get_stat_table(GopalakTNet);
end

function Result = get_stat_table(networks)
global measures
[stats,names] = get_statistics(networks);
Result = array2table(stats, 'VariableNames', measures);
Result.LeadName = names;

function [stats,names] = get_statistics(networks)
fields = fieldnames(networks);
N = length(fields);
stats = zeros(N,6);
names = cell(N,1);
for i = 1:length(fields)
    stats(i,:) = networks.(fields{i}).stat;
    names{i} = fields{i};
end