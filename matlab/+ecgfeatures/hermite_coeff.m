function [C1,C2] = hermite_coeff(Beats, I, J, L)
%
n = size(Beats,2);
C1 = zeros(n,6);
C2 = zeros(n,6);
