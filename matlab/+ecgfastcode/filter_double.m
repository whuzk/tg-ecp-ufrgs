function varargout = filter_double(varargin)

[varargout{1:nargout}] = ecgfastcode.c_filter_double(varargin{:});
clear '+ecgfastcode/c_filter_double.mexw64';