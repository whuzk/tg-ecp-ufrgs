function varargout = detect_qrs(varargin)

try
    [varargout{1:nargout}] = mex.c_detect_qrs(varargin{:});
    clear '+mex/c_detect_qrs.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_detect_qrs.mexw64';
end