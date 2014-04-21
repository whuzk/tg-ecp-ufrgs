function varargout = detect_qrs(varargin)

varargout{:} = ecgfastcode.c_detect_qrs(varargin{:});
clear '+ecgfastcode/c_detect_qrs.mexw64';