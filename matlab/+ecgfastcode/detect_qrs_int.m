function varargout = detect_qrs_int(varargin)

[varargout{1:nargout}] = ecgfastcode.c_detect_qrs_int(varargin{:});
clear '+ecgfastcode/c_detect_qrs_int.mexw64';