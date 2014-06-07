function varargout = preprocess_filter(varargin)

try
    [varargout{1:nargout}] = mexinst.c_preprocess_filter(varargin{:});
    clear '+mexinst/c_preprocess_filter.mexw64';
catch err
    disp(getReport(err));
    clear '+mexinst/c_preprocess_filter.mexw64';
end