function Result = ecg_ischemic_episode_detection(IschemicBeats, Ws)
%   Detecta episodios isquemicos com base nas batidas classificadas como
%   isquemicas e numa janela deslizante de largura especificada.
%
% Entradas:
%   IschemicBeats - batidas classificadas como isquemicas
%   Ws            - largura da janela deslizante
%
% Saída:
%   indice das batidas que determinam inicio e fim de um episodio isquemico
%
m = length(IschemicBeats);
starting = false(m,1);
ending = false(m,1);

i = Ws;
while i <= m
    window = i-Ws+1:i;
    if length(find(IschemicBeats(window))) > Ws/2
        starting(window(1)) = true;
        ending(window(end)) = true;
        i = i + Ws;
    else
        i = i + 1;
    end
end
Result = [find(starting) find(ending)];

D = Result(2:end,1) - Result(1:end-1,2);
a = find(D < Ws);
Result(a,2) = Result(a+1,2);
Result(a+1,:) = [];
