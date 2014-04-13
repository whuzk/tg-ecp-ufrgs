function Result = mmo_open_closing(Signal, B1, B2)
% Aplica a operaçao de abertura-fechamento no sinal de acordo com um par
% de elementos estruturadores B1 e B2
import ecgmmo.*;

Result = mmo_closing(mmo_opening(Signal,B1,B1),B2,B2);