function varargout = preprocess_filter(varargin)

try
    [varargout{1:nargout}] = mex.c_preprocess_filter(varargin{:});
    clear '+mex/c_preprocess_filter.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_preprocess_filter.mexw64';
end