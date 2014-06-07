function varargout = mohebbi_features(varargin)

try
    [varargout{1:nargout}] = mexinst.c_mohebbi_features(varargin{:});
    clear '+mexinst/c_mohebbi_features.mexw64';
catch err
    disp(getReport(err));
    clear '+mexinst/c_mohebbi_features.mexw64';
end