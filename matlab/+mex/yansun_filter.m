function varargout = yansun_filter(varargin)

[varargout{1:nargout}] = mex.c_yansun_filter(varargin{:});
clear '+mex/c_yansun_filter.mexw64';