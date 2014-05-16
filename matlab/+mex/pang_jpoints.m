function varargout = pang_jpoints(varargin)

try
    [varargout{1:nargout}] = mex.c_pang_jpoints(varargin{:});
    clear '+mex/c_pang_jpoints.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_pang_jpoints.mexw64';
end