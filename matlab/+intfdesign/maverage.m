function [b,a,g,d,info] = maverage(spec,varargin)
import intfdesign.*

allowed_specs = {'Width'};
specs = check_args(spec,varargin,allowed_specs);
designObj = process_specs(specs,varargin);

if isfield(designObj,'Fs')
    m = round(designObj.L * designObj.Fs);
else
    m = designObj.L;
end

if m == 0
    error('could not design a moving-average for this specfication');
end

[b,a,g,d] = design_basic_lp(1,m);

Wsl = get_sidelobe_freq(m);
info.AttenuationAtSideLobe = 20*log10(get_h(b,a,Wsl)/g);