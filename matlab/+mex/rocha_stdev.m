function varargout = rocha_stdev(varargin)

try
    [varargout{1:nargout}] = mex.c_rocha_stdev(varargin{:});
    clear '+mex/c_rocha_stdev.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_rocha_stdev.mexw64';
end