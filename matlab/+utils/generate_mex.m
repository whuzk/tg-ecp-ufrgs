blaslib = fullfile(matlabroot, 'extern', 'lib', computer('arch'), 'microsoft', 'libmwblas.lib');

%% codigos originais
mex 'ccode\c_preprocess_filter.c';
mex 'ccode\c_preprocess_detect.c';
mex 'ccode\c_mohebbi_features.c';
mex('-largeArrayDims', 'ccode\c_rocha_features.c', blaslib);
mex('-largeArrayDims', 'ccode\c_gopalak_features.c', blaslib);

%% codigos com instrumentaçao
mex 'ccode\c_preprocess_filter_inst.c';
mex 'ccode\c_preprocess_detect_inst.c';
mex 'ccode\c_mohebbi_features_inst.c';
mex('-largeArrayDims', 'ccode\c_rocha_features_inst.c', blaslib);
mex('-largeArrayDims', 'ccode\c_gopalak_features_inst.c', blaslib);

%% codigos de tempo-real
mex 'realtime\source\msfcn_baselinefit.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_beatdetector.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_beatframer.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_blockfilter.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_doublertfilter.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_edgedetector.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_intrtfilter.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_peakdetector.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_peakevaluator.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_qrsevaluator.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_rminmax.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_rrmeasurer.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_rttime.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_searchback.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_searchbest.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_searchfirst.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_searchpeak.cpp' '-Irealtime\include';
mex 'realtime\source\msfcn_segmentextractor.cpp' '-Irealtime\include';
mex('-largeArrayDims', 'realtime\source\msfcn_templatebuilder.cpp', ...
    '-Irealtime\include', blaslib);
