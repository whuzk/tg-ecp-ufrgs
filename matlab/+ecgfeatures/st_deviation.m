function [B1,B2] = st_deviation(Beats, Fs, RR)
%
n = size(Beats,2);
B1 = zeros(n,1);
B2 = zeros(n,1);
