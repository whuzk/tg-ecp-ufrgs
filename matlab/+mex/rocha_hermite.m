function varargout = rocha_hermite(varargin)

try
    [varargout{1:nargout}] = mex.c_rocha_hermite(varargin{:});
    clear '+mex/c_rocha_hermite.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_rocha_hermite.mexw64';
end