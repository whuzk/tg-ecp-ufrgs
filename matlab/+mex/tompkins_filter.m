function varargout = tompkins_filter(varargin)

try
    [varargout{1:nargout}] = mex.c_tompkins_filter(varargin{:});
    clear '+mex/c_tompkins_filter.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_tompkins_filter.mexw64';
end