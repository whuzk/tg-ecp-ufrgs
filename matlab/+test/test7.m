import nnetwork.generate_ind_stats
global basedir methods measures
path = [basedir 'extracted\'];

tic
tables = cell(8,1);
count = 0;
for i = 0:1
    disp(['Training networks for feature ' num2str(i)]);
    for j = 0:1
        disp(['Training networks for channel ' num2str(j)]);
        for k = 1:length(methods)
            if i > 0 && k < 3
                continue;
            end
            disp(['Training networks for method ' methods{k}]);
            [stats,names] = generate_ind_stats(path, i, j, methods{k}, 3);
            mytable = array2table(stats, 'VariableNames', measures);
            mytable.RecordName = names;
            count = count + 1;
            tables{count} = mytable;
        end
    end
end
toc;
save('matfiles\indstats.mat', 'tables');