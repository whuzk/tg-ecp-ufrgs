function Result = mmo_close_opening(Signal, B1, B2)
% Aplica a operaçao de fechamento-abertura no sinal de acordo com um par
% de elementos estruturadores B1 e B2
import ecgmath.*;

Result = mmo_opening(mmo_closing(Signal,B1,B1),B2,B2);