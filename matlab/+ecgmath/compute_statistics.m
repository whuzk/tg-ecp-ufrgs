function Result = compute_statistics(Known, Predicted)
% Computa as estatisticas de um diagnostico predito de ECG com base num
% diagnostico conhecido

N = length(Predicted);
if N ~= length(Known)
    error('parameter size mismatch');
end

VP = length(find( Known &  Predicted));     % verdadeiros positivos
VN = length(find(~Known & ~Predicted));     % verdadeiros negativos
FP = length(find(~Known &  Predicted));     % falsos positivos
FN = length(find( Known & ~Predicted));     % falsos negativos

Result = [
    VP/(VP+FN)      % sensibilidade
    VN/(VN+FP)      % especificidade
    VP/(VP+FP)      % preditividade positiva
    VN/(VN+FN)      % preditividade negativa
    (VP+VN)/N       % acuracia
];

% corrige casos especiais
Result(isnan(Result)) = 1;
