function varargout = mohebbi_features(varargin)

try
    [varargout{1:nargout}] = mex.c_mohebbi_features(varargin{:});
    clear '+mex/c_mohebbi_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_mohebbi_features.mexw64';
end