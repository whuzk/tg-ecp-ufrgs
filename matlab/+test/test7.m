load('matfiles\indstats.mat');
for i = 1:length(tables)
    writetable(tables{i}, ['indtable' num2str(i) '.csv']);
end
clear tables;

load('matfiles\colstats.mat');
for i = 1:length(tables)
    writetable(tables{i}, ['coltable' num2str(i) '.csv']);
end
clear tables;