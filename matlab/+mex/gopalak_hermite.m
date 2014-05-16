function varargout = gopalak_hermite(varargin)

try
    [varargout{1:nargout}] = mex.c_gopalak_hermite(varargin{:});
    clear '+mex/c_gopalak_hermite.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_gopalak_hermite.mexw64';
end