function varargout = filter_int(varargin)

[varargout{1:nargout}] = ecgfastcode.c_filter_int(varargin{:});
clear '+ecgfastcode/c_filter_int.mexw64';