function varargout = filter_int(varargin)

varargout{:} = ecgfastcode.c_filter_int(varargin{:});
clear '+ecgfastcode/c_filter_int.mexw64';