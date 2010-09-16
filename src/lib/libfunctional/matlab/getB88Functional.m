function data = getB88Functional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names{1} = 'd';
data.param_vals(1) = 0.0042;

d = sym('d');
data.functional = -d*rho_a^(4.0/3.0)*gamma_aa/(1.0+6.0*d*sqrt(gamma_aa)*asinh(sqrt(gamma_aa)))... 
                  -d*rho_b^(4.0/3.0)*gamma_bb/(1.0+6.0*d*sqrt(gamma_bb)*asinh(sqrt(gamma_bb)));
data.functional_a0 = -d*rho_b^(4.0/3.0)*gamma_bb/(1.0+6.0*d*sqrt(gamma_bb)*asinh(sqrt(gamma_bb)));
data.functional_b0 = -d*rho_a^(4.0/3.0)*gamma_aa/(1.0+6.0*d*sqrt(gamma_aa)*asinh(sqrt(gamma_aa)));
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.type = 'x';

data.name = 'B88';
data.citation = 'A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988';
data.description = 'Becke 88 Exchange (GGA Only)';
