function varargout = preprocess_detect(varargin)

try
    [varargout{1:nargout}] = mex.c_preprocess_detect(varargin{:});
    clear '+mex/c_preprocess_detect.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_preprocess_detect.mexw64';
end