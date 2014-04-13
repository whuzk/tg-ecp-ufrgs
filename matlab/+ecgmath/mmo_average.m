function Result = mmo_average(Signal, B1, B2)
% Aplica a operaçao de media no sinal de acordo com um par de elementos
% estruturadores B1 e B2
import ecgmath.*;

oc = mmo_open_closing(Signal,B1,B2);
co = mmo_close_opening(Signal,B1,B2);
Result = (oc+co)/2;