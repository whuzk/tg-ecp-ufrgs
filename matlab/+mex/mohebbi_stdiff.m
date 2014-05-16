function varargout = mohebbi_stdiff(varargin)

try
    [varargout{1:nargout}] = mex.c_mohebbi_stdiff(varargin{:});
    clear '+mex/c_mohebbi_stdiff.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_mohebbi_stdiff.mexw64';
end