blaslib = fullfile(matlabroot, 'extern', 'lib', computer('arch'), 'microsoft', 'libmwblas.lib');

mex 'ccode\c_preprocess_filter.c';
mex 'ccode\c_preprocess_detect.c';
mex 'ccode\c_pang_jpoints.c';
mex 'ccode\c_mohebbi_ijpoints.c';
mex 'ccode\c_mohebbi_stdiff.c';
mex 'ccode\c_mohebbi_segments.c';
mex 'ccode\c_rocha_stdev.c';
mex 'ccode\c_rocha_segments.c';
mex 'ccode\c_gopalak_segments.c';
mex('-largeArrayDims', 'ccode\c_rocha_hermite.c', blaslib);
mex('-largeArrayDims', 'ccode\c_gopalak_hermite.c', blaslib);