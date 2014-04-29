function R = adjust_qrs(data,lead,R,Fs)

if ismember(lead,{'MLIII'})
    data = -data;
end

L = floor(0.05*Fs);
R = ecgmath.neighbour_max(data,R,L);