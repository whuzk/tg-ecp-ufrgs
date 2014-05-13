function Result = get_m(N,Wc,sense)

switch sense
    case '3db'
        C = 1/sqrt(2);
        Result = (0.91823/sqrt(N+0.072177))/Wc;
    case '6db'
        C = 1/2;
        Result = (1.2992/sqrt(N+0.15023))/Wc;
    case '24db'
        C = 1-1/sqrt(2);
        Result = (1.732/sqrt(N+0.28356))/Wc;
    case 'nom'
        Result = 2/Wc;
end
%{
if ~strcmp(sense,'nom')
    syms m;
    w = Wc*pi/2;
    eq = sym((sin(m*w)/sin(w)/m)^N == C);
    Result = eval(vpasolve(eq, m, Result));
end
%}