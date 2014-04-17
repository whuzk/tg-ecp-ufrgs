function Result = compute_statistics(Known, Predicted)
% Computa as estatisticas de um diagnostico predito de ECG com base num
% diagnostico conhecido

if length(Known) ~= length(Predicted)
    error('parameter size mismatch');
end

VP = length(find( Known &  Predicted));     % verdadeiros positivos
VN = length(find(~Known & ~Predicted));     % verdadeiros negativos
FP = length(find(~Known &  Predicted));     % falsos positivos
FN = length(find( Known & ~Predicted));     % falsos negativos
TT = length(find(Known));                   % total conhecido

Result = [
    VP/(VP+FN)      % sensibilidade
    VN/(VN+FP)      % especificidade
    VP/(VP+FP)      % preditividade positiva
    VN/(VN+FN)      % preditividade negativa
    (VP+VN)/TT      % acuracia
    (FP+FN)/TT      % taxa de detecçcao falsa
];