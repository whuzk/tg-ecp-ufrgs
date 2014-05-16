function varargout = extract_beats(varargin)

try
    [varargout{1:nargout}] = mex.c_extract_beats(varargin{:});
    clear '+mex/c_extract_beats.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_extract_beats.mexw64';
end