function R = adjust_qrs(data,lead,R,Fs)
% Funçao para ajustar a localizaçao dos picos de onda R, de acordo com a
% derivaçao
global rwavepolarity

data = rwavepolarity(lead) * data;
R = utils.neighbour_max(data,R,floor(0.05*Fs));