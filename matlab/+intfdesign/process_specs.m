function designObj = process_specs(specs, args)

if length(args) > length(specs)
    designObj.Fs = args{end};
end

for i = 1:length(specs)
    switch specs{i}
        case 'N'
            designObj.N = args{i};
        case 'Width'
            designObj.L = args{i};
        case 'M'
            designObj.M = args{i};
        case {'F3db','F6db','Fnom'}
            designObj.Fc = args{i};
            designObj.sense = strrep(specs{i},'F','');
        case {'BW3db','BW6db','BWnom'}
            designObj.Fc = args{i}(1);
            designObj.Bw = args{i}(2);
            designObj.sense = strrep(specs{i},'BW','');
    end
end