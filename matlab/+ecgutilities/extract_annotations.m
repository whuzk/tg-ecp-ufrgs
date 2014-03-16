function Result = extract_annotations(text)

if strfind(text, '[')
    annot = textscan(text, '%s%s%s%s%s%s%s%s');
    annot{1} = strcat(annot{1}, {' '}, annot{2});
    annot(2) = [];
elseif strfind(text, 'Time')
    annot = textscan(text, '%s%s%s%s%s%s%s', 'HeaderLines', 1);
else
    annot = textscan(text, '%s%s%s%s%s%s%s');
end

Result.DateTime = annot{1};
Result.Sample = annot{2};
Result.Type = annot{3};
Result.Sub = annot{4};
Result.Chan = annot{5};
Result.Num = annot{6};
Result.Aux = annot{7};
