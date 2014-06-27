function plot_class_distribution(count)
global stclasses tclasses classes

fields = fieldnames(count);
for i = 1:length(fields)
    C = count.(fields{i});
    C1 = reshape(C,1,numel(C));
    C2 = sum(C,1);
    C3 = sum(C,2);
    figure, plot_pie(C1, classes); title(fields{i});
    figure, plot_pie(C2, tclasses); title(fields{i});
    figure, plot_pie(C3, stclasses); title(fields{i});
    pause;
end

function plot_pie(C,classes)
total = sum(C);
idx = C./total < 0.02;
if length(find(idx)) > 1
    C = [C(~idx) sum(C(idx))];
    L = {classes{~idx} 'other'};
else
    L = classes;
end
pie(C(C > 0));
legend(L{C > 0},'Location','southeast');