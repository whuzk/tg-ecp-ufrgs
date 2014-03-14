Records = fieldnames(Beatset);
k = numel(Records);
temp = zeros(k,1);
for i = 1:k
    Data = Beatset.(Records{i});
    temp(i) = length(Data.R);
end
temp