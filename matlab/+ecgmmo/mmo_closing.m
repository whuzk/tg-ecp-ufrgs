function Result = mmo_closing(Signal, B1, B2)
% Aplica a operaçao de fechamento no sinal de acordo com um par de
% elementos estruturadores B1 e B2
import ecgmmo.*;

Result = mmo_erosion(mmo_dilation(Signal,B1),B2);