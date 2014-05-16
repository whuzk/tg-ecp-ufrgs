function varargout = rocha_segments(varargin)

try
    [varargout{1:nargout}] = mex.c_rocha_segments(varargin{:});
    clear '+mex/c_rocha_segments.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_rocha_segments.mexw64';
end