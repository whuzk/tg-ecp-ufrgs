function [b,a,g,d,info] = lowpass(spec,varargin)
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

if Wc == 0
    error('could not design a low-pass for this specfication');
elseif Wc == 1
    % design the simplest low-pass
    [b,a,g,d] = design_basic_lp(N,2);
    Wsl = -Inf;
elseif Wc <= 0.5 || strcmp(sense,'nom')
    % design low-pass directly
    m = ceil(get_m(N,Wc,sense));
    [b,a,g,d] = design_basic_lp(N,m);
    Wsl = get_sidelobe_freq(m);
else
    % subtract high-pass from all-pass
    sense = reciprocal_sense(sense);
    m = floor(get_m(N,1-Wc,sense));
    if mod(floor(m/2),2) ~= 0
        error('could not design a low-pass for this specfication');
    end
    [b,a,g,d] = design_basic_hp(N,m);
    [b,d] = subtract_allpass(b,a,g,d);
    Wsl = -Inf;
end

info.AttenuationAtDesiredFc = 20*log10(get_h(b,a,Wc)/g);
info.AttenuationAtSideLobe = 20*log10(get_h(b,a,Wsl)/g);