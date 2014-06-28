load('matfiles\indstats.mat');
lists = cell(length(tables),1);
for i = 1:length(tables)
    mytable = tables{i};
    sel1 = mytable.Sensitivity >= 0.8;
    sel2 = mytable.PositivePred >= 0.8;
    lists{i} = mytable(sel1 & sel2, 'RecordName');
end
save('matfiles\selection.mat', 'lists');

merge = cell(4,1);
for i = 1:2
    list1 = lists{3*(i-1)+1};
    list2 = lists{3*(i-1)+2};
    list3 = lists{3*(i-1)+3};
    list4 = union(union(list1{:,:},list2{:,:}),list3{:,:});
    idx = true(size(list4));
    for j = 1:length(list4)
        name = list4{j};
        cond1 = any(strcmp(name, list1{:,:}));
        cond2 = any(strcmp(name, list2{:,:}));
        cond3 = any(strcmp(name, list3{:,:}));
        if cond1 + cond2 + cond3 <= 1
            idx(j) = false;
        end
    end
    merge{i} = array2table(list4(idx),'VariableNames',{'RecordName'});
end
merge{3} = lists{end-1};
merge{4} = lists{end};
save('matfiles\merge.mat', 'merge');