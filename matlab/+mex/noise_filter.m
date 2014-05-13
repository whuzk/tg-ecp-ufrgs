function varargout = noise_filter(varargin)

try
    [varargout{1:nargout}] = mex.c_noise_filter(varargin{:});
    clear '+mex/c_noise_filter.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_noise_filter.mexw64';
end