function varargout = detect_qrs(varargin)

[varargout{1:nargout}] = ecgfastcode.c_detect_qrs(varargin{:});
clear '+ecgfastcode/c_detect_qrs.mexw64';