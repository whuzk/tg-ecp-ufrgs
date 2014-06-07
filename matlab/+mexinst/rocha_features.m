function varargout = rocha_features(varargin)

try
    [varargout{1:nargout}] = mexinst.c_rocha_features(varargin{:});
    clear '+mexinst/c_rocha_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mexinst/c_rocha_features.mexw64';
end