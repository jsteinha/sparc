function sig=getsig(s1,s2,m)
    global u0 s0 dof dof2;
    if m > 0
        sig = s0 + s2 - s1*s1'/m + dof2*m/(m+dof2)*(s1/m-u0)*(s1/m-u0)';
    else
        sig = s0;
    end
    sig = sig/(m+dof);
end