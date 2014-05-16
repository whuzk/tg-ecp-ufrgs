function varargout = detect_qrs_prod(varargin)

try
    [varargout{1:nargout}] = mex.c_detect_qrs_prod(varargin{:});
    clear '+mex/c_detect_qrs_prod.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_detect_qrs_prod.mexw64';
end