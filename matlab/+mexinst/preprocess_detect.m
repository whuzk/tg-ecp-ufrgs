function varargout = preprocess_detect(varargin)

try
    [varargout{1:nargout}] = mexinst.c_preprocess_detect(varargin{:});
    clear '+mexinst/c_preprocess_detect.mexw64';
catch err
    disp(getReport(err));
    clear '+mexinst/c_preprocess_detect.mexw64';
end