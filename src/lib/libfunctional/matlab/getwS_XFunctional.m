function data = getwS_XFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names{1} = 'omega';
data.param_vals(1) = 0.3; 

syms omega real

Ka = 3*(3/(4*pi))^(1/3);
Kb = 3*(3/(4*pi))^(1/3);
aa = omega * Ka^(1/2) / (6 * sqrt(pi) * rho_a^(1/3));
ab = omega * Kb^(1/2) / (6 * sqrt(pi) * rho_b^(1/3));

ba = exp(-1/(4 * aa^2)) - 1;
bb = exp(-1/(4 * ab^2)) - 1;

ca = 2*aa^2 * ba - 1;
cb = 2*ab^2 * bb - 1;

Aa = 1 - 8/3*aa*(sqrt(pi) * erf(1/(2*aa)) + 2*aa*(ba - ca));
Ab = 1 - 8/3*ab*(sqrt(pi) * erf(1/(2*ab)) + 2*ab*(bb - cb));

data.functional = - 1/2 * (Ka * Aa * rho_a^(4.0/3.0) + Kb * Ab * rho_b^(4.0/3.0)); 
data.functional_a0 = - 1/2 * (Ka * Ab * rho_b^(4.0/3.0)); 
data.functional_b0 = - 1/2 * (Kb * Aa * rho_a^(4.0/3.0)); 
data.functional_a0b0 = 0; 

data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 1;

data.name = 'wS_X';
data.citation = 'Null';
data.description = 'Range-Corrected Slater LSDA';
