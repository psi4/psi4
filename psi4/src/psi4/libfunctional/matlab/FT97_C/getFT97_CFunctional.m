function data = getFT97_CFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c0', 'c', 'tspi_m13',...
    'a1', 'a2', 'a3', 'a4', 'a5',...
    'kaa0', 'kaa1', 'kaa2', 'raa1', 'raa2',...
    'kab0', 'kab1', 'rab1', 'k1', 'k2'};

data.param_vals = [(1.0-log(2.0))/(pi*pi), 1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), (36*pi)^(-1.0/3.0)...
    1.622118767, 0.489958076, 1.379021941, 4.946281353, 3.600612059,...
    1.200801774,-0.812904345, 0.859614445, 0.655638823, 1.089338848,...
    1.291551074,-0.349064173, 0.083275880, 0.939016000, 1.733170000];
syms c0 c tspi_m13 real
syms a1 a2 a3 a4 a5 real
syms kaa0 kaa1 kaa2 raa1 raa2 real
syms kab0 kab1 rab1 k1 k2 real
syms r chi x real

Eii = int(exp(x)/x,x);
rs = c*r^(-1/3); 
chi_a = sqrt(gamma_aa)/rho_a^(4/3);
chi_b = sqrt(gamma_bb)/rho_b^(4/3);

Fab = (1 + a1 * (-tspi_m13*chi)^2 + a2^2 * (-tspi_m13 * chi)^4) * exp (-a2^2 * (-tspi_m13*chi)^4) / ...
    sqrt(1 + a3* ( -tspi_m13 * chi)^2 / rs); 
Faa = (1 + a4^2 * (-tspi_m13 * chi)^4) * exp (-a4^2 * (-tspi_m13 * chi)^4) / ...
    sqrt(1 + a5* ( -tspi_m13 * chi)^2 / rs); 

kaa = kaa0 + kaa1 * (1 - exp(-raa1*rs^(2/5))) + kaa2 * (1- exp(-raa2*rs^(1/2)));
kab = kab0 + kab1 * (1 - exp(-rab1*rs^(4/5)));

rc_aa = kaa*Faa;
rc_ab = kab*Fab;

uab = 2/3*c0*rs/(rc_ab)^2; 
uaa = 2/3*c0*rs/(rc_aa)^2;

Aab = (6 + 4*sqrt(uab) + 4*uab) / (3 + 6*sqrt(uab) + 6*uab); 
Aaa = (6 + 4*sqrt(uaa) + 4*uaa) / (3 + 6*sqrt(uaa) + 6*uaa);

F = exp(-rs^2/(k1*sqrt(rs) + k2*rs)^2);


eab = 1/2*c0 * (exp(uab)*subs(Eii,x,-uab) + Aab*(uab*exp(uab)*subs(Eii,x,-uab)+1)) * heaviside(100000*rc_ab - 2/3*c0*rs);  
eaa = 1/2*F*c0 * (exp(uaa)*subs(Eii,x,-uaa) + Aaa*(uaa*exp(uaa)*subs(Eii,x,-uaa)+1)) * heaviside(100000*rc_ab - 2/3*c0*rs); 

data.functional = rho_a * (subs(eaa,{r,chi}, {rho_a, chi_a}) + subs(eab, {r, chi}, {rho_b, chi_b})) + ... 
                   rho_b * (subs(eaa,{r,chi}, {rho_b, chi_b}) + subs(eab, {r, chi}, {rho_a, chi_a}));
data.functional_a0 = rho_b*subs(eaa,{r,chi},{rho_b,chi_b}); 
data.functional_b0 = rho_a*subs(eaa,{r,chi},{rho_a,chi_a}); 
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'FT97_C';
data.citation = 'M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997';
data.description = 'FT97 Correlation (Involves Ei functions)'; 
