function data = getS_XFunctional()

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

c = sym('c');
data.functional = - c*(rho_a^(4.0/3.0)+rho_b^(4.0/3.0)); 
data.functional_a0 = - c*(rho_b^(4.0/3.0)); 
data.functional_b0 = - c*(rho_a^(4.0/3.0)); 
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 1;

data.name = 'S_X';
data.citation = 'J.C. Slater, Phys. Rev., 81(3):385-390, 1951';
data.description = 'Simple Slater LSDA exchange functional';
