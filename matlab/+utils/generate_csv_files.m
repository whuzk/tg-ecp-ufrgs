function generate_csv_files(prefix,tables)

avgtable = table;
for i = 1:length(tables)
    newtable = average_table(['Table' num2str(i)], tables{i});
    avgtable = [avgtable; newtable];
end
firstheader = tables{1}.Properties.VariableNames(1);
lastrow = table({'Average'}, 'VariableNames', firstheader);
firsttable = [tables{1}; [lastrow avgtable(1,2:end)]];
writetable(format_table(firsttable), [prefix 'table1.csv']);
writetable(format_table(avgtable), [prefix 'average.csv']);


function outtable = format_table(intable)
outtable = table;
column = intable.Properties.VariableNames{1};
outtable.(column) = intable{:,1};
for i = 2:7
    column = intable.Properties.VariableNames{i};
    outtable.(column) = num2str(intable{:,i}*100,'%.1f');
end

function outtable = average_table(name,intable)
outtable = table;
outtable.TableName = {name};
for i = 2:7
    column = intable.Properties.VariableNames{i};
    outtable.(column) = nanmean(intable{:,i});
end