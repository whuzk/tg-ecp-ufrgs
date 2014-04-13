function Result = mmo_opening(Signal, B1, B2)
% Aplica a operaçao de abertura no sinal de acordo com um par de elementos
% estruturadores B1 e B2
import ecgmmo.*;

Result = mmo_dilation(mmo_erosion(Signal,B1),B2);