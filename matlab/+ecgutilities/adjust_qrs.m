function R = adjust_qrs(data,lead,R,Fs)
% Fun�ao para ajustar a localiza�ao dos picos de onda R, de acordo com a
% deriva�ao
global rwavepolarity

data = rwavepolarity(lead) * data;
R = utils.neighbour_max(data,R,floor(0.05*Fs));