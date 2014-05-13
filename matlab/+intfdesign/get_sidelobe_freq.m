function Result = get_sidelobe_freq(m)

if m > 2
    syms u;
    a = m+1;
    b = m-1;
    eq = sym(b*sin(a*u*pi/2) == a*sin(b*u*pi/2));
    Result = eval(vpasolve(eq, u, 2.8607/m));
else
    Result = -Inf;
end