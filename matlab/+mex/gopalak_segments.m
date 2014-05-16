function varargout = gopalak_segments(varargin)

try
    [varargout{1:nargout}] = mex.c_gopalak_segments(varargin{:});
    clear '+mex/c_gopalak_segments.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_gopalak_segments.mexw64';
end