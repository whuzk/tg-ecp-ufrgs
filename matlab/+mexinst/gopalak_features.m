function varargout = gopalak_features(varargin)

try
    [varargout{1:nargout}] = mexinst.c_gopalak_features(varargin{:});
    clear '+mexinst/c_gopalak_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mexinst/c_gopalak_features.mexw64';
end