function Result = sun_filter(Signal,Fs)

Temp1 = ecgfilter.krishnan_filter(Signal,Fs);
Result = ecgmmo.mmo_derivative_flat(Temp1,20);