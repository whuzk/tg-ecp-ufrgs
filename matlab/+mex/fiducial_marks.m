function Result = fiducial_marks(varargin)

[A,B,C,D,E,F,G,I,J] = mex.c_fiducial_marks(varargin{:});
clear '+mex/c_fiducial_marks.mexw64';

Result = struct('P',[A B],'R',[C D E],'T',[F G],'IJ',[I J]);