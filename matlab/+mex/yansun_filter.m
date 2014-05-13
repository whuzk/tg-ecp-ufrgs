function varargout = yansun_filter(varargin)

try
    [varargout{1:nargout}] = mex.c_yansun_filter(varargin{:});
    clear '+mex/c_yansun_filter.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_yansun_filter.mexw64';
end