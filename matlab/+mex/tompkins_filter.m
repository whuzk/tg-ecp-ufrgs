function varargout = tompkins_filter(varargin)

[varargout{1:nargout}] = mex.c_tompkins_filter(varargin{:});
clear '+mex/c_tompkins_filter.mexw64';