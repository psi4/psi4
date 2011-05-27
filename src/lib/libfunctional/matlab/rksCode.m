function rksCode(functional, name)

syms rho_a rho_b gamma_aa gamma_ab gamma_bb tau_a tau_b real;
functional = subs(functional,rho_b,rho_a,0);
functional = subs(functional,gamma_ab,gamma_aa,0);
functional = subs(functional,gamma_bb,gamma_aa,0);
functional = subs(functional,tau_b,tau_a,0);

try
    syms qqqq;
    ccode(subs(functional + qqqq, qqqq, 0),'file',name);
catch exception
    syms qqqq;
    ccode(diff(functional, qqqq),'file',name);
end
