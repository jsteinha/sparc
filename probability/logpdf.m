function E = logpdf(x,s1,s2,m)
    global conc dof;
    %mu = (s1+dof2*u0)/(m+dof2);
    %if m > 0
    %    sig = s0 + s2 - s1*s1'/m + dof2*m/(m+dof2)*(s1/m-u0)*(s1/m-u0)';
    %else
    %    sig = s0;
    %end
    %sig = sig/(m+dof);
    mu = getmu(s1,m);
    sig = getsig(s1,s2,m);
    E = log(mvtpdf(x-mu,sig,m+dof))+log(m+conc);
end
