function varargout = detect_qrs_double(varargin)

[varargout{1:nargout}] = ecgfastcode.c_detect_qrs_double(varargin{:});
clear '+ecgfastcode/c_detect_qrs_double.mexw64';