function varargout = template_and_artifacts(varargin)

try
    [varargout{1:nargout}] = mex.c_template_and_artifacts(varargin{:});
    clear '+mex/c_template_and_artifacts.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_template_and_artifacts.mexw64';
end