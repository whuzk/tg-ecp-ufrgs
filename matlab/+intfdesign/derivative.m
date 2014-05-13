function [b,a,g,d,info] = derivative(spec,varargin)
import intfdesign.*

allowed_specs = {'N,M'};
specs = check_args(spec,varargin,allowed_specs);
designObj = process_specs(specs,varargin);

N = designObj.N;
if N == 0
    error('N must be positive for design of a derivative filter');
end
M = designObj.M;

[b,a,g,d] = design_basic_de(N,M);

info = [];