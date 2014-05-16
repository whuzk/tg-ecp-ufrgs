function varargout = mohebbi_ijpoints(varargin)

try
    [varargout{1:nargout}] = mex.c_mohebbi_ijpoints(varargin{:});
    clear '+mex/c_mohebbi_ijpoints.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_mohebbi_ijpoints.mexw64';
end