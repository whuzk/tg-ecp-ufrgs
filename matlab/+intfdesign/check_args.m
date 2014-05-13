function specs = check_args(spec,args,allowed)

if ~ismember(spec,allowed)
    error('specification string ''%s'' unknown for this filter type', spec);
end

specs = strsplit(spec, ',');

N = length(args);
M = length(specs);
if N < M || N > M+1
    error('wrong number of arguments for the given specification string');
end