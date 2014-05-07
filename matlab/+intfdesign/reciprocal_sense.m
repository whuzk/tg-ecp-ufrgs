function Result = reciprocal_sense(sense)

if strcmp(sense,'3db')
    Result = '24db';
elseif strcmp(sense,'24db')
    Result = '3db';
else
    Result = sense;
end