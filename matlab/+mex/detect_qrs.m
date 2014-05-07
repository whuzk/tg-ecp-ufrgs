function varargout = detect_qrs(varargin)

[varargout{1:nargout}] = mex.c_detect_qrs(varargin{:});
clear '+mex/c_detect_qrs.mexw64';