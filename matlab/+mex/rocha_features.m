function varargout = rocha_features(varargin)

try
    [varargout{1:nargout}] = mex.c_rocha_features(varargin{:});
    clear '+mex/c_rocha_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_rocha_features.mexw64';
end