function Result = fiducial_marks(varargin)

try
    [A,B,C,D,E,F,G] = mex.c_fiducial_marks(varargin{:});
    Result = [A(:) B(:) C(:) D(:) E(:) F(:) G(:)]';
    clear '+mex/c_fiducial_marks.mexw64';
catch err
    disp(getReport(err));
    clear '+mex/c_fiducial_marks.mexw64';
end