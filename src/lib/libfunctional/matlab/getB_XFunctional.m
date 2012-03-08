function data = getB_XFunctional()

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

syms d real

K = 3*(3/(4*pi))^(1/3);

chi_a = sqrt(gamma_aa) / rho_a^(4.0/3.0); 
chi_b = sqrt(gamma_bb) / rho_b^(4.0/3.0); 

Ka = (K + 2 * d  * chi_a^2 / (1 + 6 * d * chi_a * asinh(chi_a)));
Kb = (K + 2 * d  * chi_b^2 / (1 + 6 * d * chi_b * asinh(chi_b)));

data.functional = - 1/2 * (Ka * rho_a^(4.0/3.0) + Kb * rho_b^(4.0/3.0)); 
data.functional_a0 = - 1/2 * (Kb * rho_b^(4.0/3.0)); 
data.functional_b0 = - 1/2 * (Ka * rho_a^(4.0/3.0)); 
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 1;
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.type = 'x';

data.name = 'B_X';
data.citation = 'A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988';
data.description = 'Becke Exchange (S+B88)';
