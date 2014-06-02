blaslib = fullfile(matlabroot, 'extern', 'lib', computer('arch'), 'microsoft', 'libmwblas.lib');

mex 'ccode\c_preprocess_filter.c';
mex 'ccode\c_preprocess_detect.c';
mex 'ccode\c_mohebbi_features.c';
mex('-largeArrayDims', 'ccode\c_rocha_features.c', blaslib);
mex('-largeArrayDims', 'ccode\c_gopalak_features.c', blaslib);