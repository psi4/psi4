function data = getB88Functional()

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
data.functional = - c*(rho_a^(4.0/3.0)+rho_b^(4.0/3.0))...
    -d*rho_a^(4.0/3.0)*gamma_aa/(1.0+6.0*d*sqrt(gamma_aa)*asinh(sqrt(gamma_aa)))... 
    -d*rho_b^(4.0/3.0)*gamma_bb/(1.0+6.0*d*sqrt(gamma_bb)*asinh(sqrt(gamma_bb)));

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 1;

data.name = 'B88';
data.citation = 'A.D. Becke, Phys. Rev. A, 38(6):3098–3100, 1988';
data.description = 'Becke 88 Exchange';
