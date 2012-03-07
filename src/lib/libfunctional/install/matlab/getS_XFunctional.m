function data = getS_XFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

Ka = 3*(3/(4*pi))^(1/3);
Kb = 3*(3/(4*pi))^(1/3);

data.functional = - 1/2 * (Ka * rho_a^(4.0/3.0) + Kb * rho_b^(4.0/3.0)); 
data.functional_a0 = - 1/2 * (Ka * rho_b^(4.0/3.0)); 
data.functional_b0 = - 1/2 * (Kb * rho_a^(4.0/3.0)); 
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 1;

data.name = 'S_X';
data.citation = 'J.C. Slater, Phys. Rev., 81(3):385-390, 1951';
data.description = 'Simple Slater LSDA exchange functional';
