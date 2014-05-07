function [b,a,g,d,info] = bandpass(spec,varargin)
import intfdesign.*

allowed_specs = {'N,F3db' 'N,F6db' 'N,BW3db' 'N,BW6db'};
specs = check_args(spec,varargin,allowed_specs);
designObj = process_specs(specs,varargin);

N = designObj.N;
if isfield(designObj,'Fs')
    Wc = designObj.Fc * 2/designObj.Fs;
else
    Wc = designObj.Fc;
end
sense = designObj.sense;

if isfield(designObj,'Bw')
    centerFreq = Wc;
    if isfield(designObj,'Fs')
        Bw = designObj.Bw * 2/designObj.Fs;
    else
        Bw = designObj.Bw;
    end
    Wc = [Wc-Bw/2 Wc+Bw/2];
else
    centerFreq = sum(Wc)/2;
    Bw = diff(Wc);
end

% try to design the bandpass bandpass
[p,q] = rat(centerFreq/2,1e-2);
if (p == 1) && (q == 3 || q == 4 || q == 6)
    m = get_m(N,Bw/2,sense);
    m = ceil(m/q)*q;
    [b,a,g,d] = design_basic_bp(N,m,pi*2/q);
    Wsl = 2/q + get_sidelobe_freq(m)*[-1 1];
    info.AttenuationAtSideLobe = 20*log10(get_h(b,a,Wsl)/g);
    info.AttenuationAtSideLobe(Wsl < 0 | Wsl > 1) = -Inf;
else
    error('could not design a bandpass for this specfication');
end

info.AttenuationAtDesiredFc = 20*log10(get_h(b,a,Wc)/g);