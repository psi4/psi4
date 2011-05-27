function data = getB_XFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names{1} = 'c';
data.param_vals(1) = 3/8*3^(1/3)*4^(2/3)*pi^(-1/3);
data.param_names{2} = 'd';
data.param_vals(2) = 0.0042;

c = sym('c');
d = sym('d');
chi_a = sym('chi_a');
chi_b = sym('chi_b');
data.functional = -c*(rho_a^(4.0/3.0)+rho_b^(4.0/3.0)) ... 
                  -d*rho_a^(4.0/3.0)*chi_a^2/(1.0+6.0*d*chi_a*asinh(chi_a))... 
                  -d*rho_b^(4.0/3.0)*chi_b^2/(1.0+6.0*d*chi_b*asinh(chi_b));
data.functional = subs(data.functional,chi_a,sqrt(gamma_aa)/rho_a^(4.0/3.0),0);
data.functional = subs(data.functional,chi_b,sqrt(gamma_bb)/rho_b^(4.0/3.0),0);
data.functional_a0 = -c*(rho_b^(4.0/3.0)) ...
                  -d*rho_b^(4.0/3.0)*chi_b^2/(1.0+6.0*d*chi_b*asinh(chi_b));
data.functional_a0 = subs(data.functional_a0,chi_b,sqrt(gamma_bb)/rho_b^(4.0/3.0),0);
data.functional_b0 = -c*(rho_a^(4.0/3.0)) ...
                  -d*rho_a^(4.0/3.0)*chi_a^2/(1.0+6.0*d*chi_a*asinh(chi_a));
data.functional_b0 = subs(data.functional_b0,chi_a,sqrt(gamma_aa)/rho_a^(4.0/3.0),0);
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.type = 'x';

data.name = 'B_X';
data.citation = 'A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988';
data.description = 'Becke Exchange (S+B88)';
