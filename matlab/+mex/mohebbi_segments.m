function varargout = mohebbi_segments(varargin)

try
    [varargout{1:nargout}] = mex.c_mohebbi_segments(varargin{:});
    clear '+mex/c_mohebbi_segments.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_mohebbi_segments.mexw64';
end