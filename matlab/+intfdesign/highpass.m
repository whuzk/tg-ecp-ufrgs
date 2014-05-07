function [b,a,g,d,info] = highpass(spec,varargin)
import intfdesign.*

allowed_specs = {'N,F3db' 'N,F6db' 'N,F24db' 'N,Fnom'};
specs = check_args(spec,varargin,allowed_specs);
designObj = process_specs(specs,varargin);

N = designObj.N;
if isfield(designObj,'Fs')
    Wc = designObj.Fc * 2/designObj.Fs;
else
    Wc = designObj.Fc;
end
sense = designObj.sense;

if Wc == 1
    error('could not design a high-pass for this specfication');
elseif Wc == 0
    % design the simplest high-pass
    [b,a,g,d] = design_basic_hp(N,2);
    Wsl = -Inf;
elseif Wc >= 0.5 || strcmp(sense,'nom')
    % design high-pass directly
    m = ceil(get_m(N,1-Wc,sense));
    [b,a,g,d] = design_basic_hp(N,m);
    Wsl = 1 - get_sidelobe_freq(m);
else
    % subtract low-pass from all-pass
    sense = reciprocal_sense(sense);
    m = floor(get_m(N,Wc,sense));
    [b,a,g,d] = design_basic_lp(N,m);
    [b,d] = subtract_allpass(b,a,g,d);
    Wsl = -Inf;
end

info.AttenuationAtDesiredFc = 20*log10(get_h(b,a,Wc)/g);
info.AttenuationAtSideLobe = 20*log10(get_h(b,a,Wsl)/g);