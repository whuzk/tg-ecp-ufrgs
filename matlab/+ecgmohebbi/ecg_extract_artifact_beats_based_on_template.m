function Result = ecg_extract_artifact_beats_based_on_template(Signal, R, Template)
%   Detecta batidas anormais presentes no sinal de ECG, com base num
%   template de batida normal.
%
% Entradas:
%   Signal   - amplitudes normalizadas do sinal
%   R        - localizaçao dos picos de onda R
%   Template - template de onda normal
%
% Saída:
%   indices das batidas classificadas como anormais
%
N = length(Signal);
m = size(R,1);
l = length(Template);
l1 = fix(l/2);
l2 = l - l1;

Result = false(m,1);
for i = 1:m
    left = max(1,R(i)-l1);
    right = min(N,R(i)+l2-1);
    Beat = Signal(left:right);
    samples = left+l1-R(i)+(1:length(Beat));
    Difference = Beat-Template(samples);
    Result(i) = norm(Difference) > 1.0;
end
