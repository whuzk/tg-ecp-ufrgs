function varargout = gopalak_features(varargin)

try
    [varargout{1:nargout}] = mex.c_gopalak_features(varargin{:});
    clear '+mex/c_gopalak_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_gopalak_features.mexw64';
end