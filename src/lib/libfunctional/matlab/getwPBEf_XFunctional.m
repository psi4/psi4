function data = getwPBEf_XFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'k','e','kp','mu_'};
data.param_vals = [3.0^(1.0/3.0)*pi^(2.0/3.0), -3.0/4.0/pi, 0.804, 0.2195149727645171];

syms k e kp mu_ real;
syms n s real;

kF = k*n^(1/3);
eX = e*kF;
S = 1/2*sqrt(s)/(kF*n);

FX = 1 + kp - kp*(1 + mu_*S^2/kp)^-1;

data.functional = rho_a*subs(eX,n,2*rho_a)*subs(FX,{n,s},{2*rho_a,4*gamma_aa}) + ... 
                  rho_b*subs(eX,n,2*rho_b)*subs(FX,{n,s},{2*rho_b,4*gamma_bb});
data.functional_a0 = rho_b*subs(eX,n,2*rho_b)*subs(FX,{n,s},{2*rho_b,4*gamma_bb}); 
data.functional_b0 = rho_a*subs(eX,n,2*rho_a)*subs(FX,{n,s},{2*rho_a,4*gamma_aa});  
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.type = 'x';

data.name = 'wPBEf_X';
data.citation = 'J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996';
data.description = 'Fake SR PBE Exchange Functional';
