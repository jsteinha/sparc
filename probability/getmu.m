function mu=getmu(s1,m)
    global u0 dof2;
    mu = (s1+dof2*u0)/(m+dof2);
end
