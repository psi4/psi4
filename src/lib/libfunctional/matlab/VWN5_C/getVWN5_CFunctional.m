function data = getVWN5_CFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names{1} = 'c';
data.param_vals(1) = 1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*(1.0/pi)^(1.0/3.0);
% d2fz0 = d_zz[fz] (0) 
data.param_names{2} = 'd2fz0';
data.param_vals(2) = 4.0/9.0/(2.0^(1.0/3.0)-1.0);

% EcP parameters
data.param_names{3} = 'EcP_1';
data.param_vals(3) = 0.03109070000;
data.param_names{4} = 'EcP_2';
data.param_vals(4) = -0.10498;
data.param_names{5} = 'EcP_3';
data.param_vals(5) = 3.72744;
data.param_names{6} = 'EcP_4';
data.param_vals(6) = 12.9352;

% EcF parameters
data.param_names{7} = 'EcF_1';
data.param_vals(7) = 0.01554535000;
data.param_names{8} = 'EcF_2';
data.param_vals(8) = -0.32500;
data.param_names{9} = 'EcF_3';
data.param_vals(9) = 7.06042;
data.param_names{10} = 'EcF_4';
data.param_vals(10) = 18.0578;

% Ac parameters
data.param_names{11} = 'Ac_1';
data.param_vals(11) = -1.0/6.0*pi^(-2.0);
data.param_names{12} = 'Ac_2';
data.param_vals(12) = -0.00475840;
data.param_names{13} = 'Ac_3';
data.param_vals(13) = 1.13107;
data.param_names{14} = 'Ac_4';
data.param_vals(14) = 13.0045;

% 2^(1/3)
data.param_names{15} = 'two_13';
data.param_vals(15) = 2.0^(1.0/3.0);

syms r c d2fz0 two_13 real
syms EcP_1 EcP_2 EcP_3 EcP_4 real
syms EcF_1 EcF_2 EcF_3 EcF_4 real
syms Ac_1 Ac_2 Ac_3 Ac_4 real

rho = rho_a + rho_b;
z = (rho_a - rho_b)/(rho);
fz = ((1+r)^(4/3) + (1-r)^(4/3) - 2)/(2*two_13 - 2);
rs = c/rho^(1/3);
x = sqrt(rs);

syms B C real
X = x^2 + B*x + C;
Q = sqrt(4*C-B^2);

syms A X0 real

E = A*(log(x^2/X) + 2*B*atan(Q/(2*x+B))/Q - B*X0*...
    (log((x-X0)^2/X) + 2*(B+2*X0)*atan(Q/(2*x+B))/Q)/(X0^2 + B*X0 + C));

EcP = subs(E,{A,X0,B,C},{EcP_1,EcP_2,EcP_3,EcP_4});
EcF = subs(E,{A,X0,B,C},{EcF_1,EcF_2,EcF_3,EcF_4});
Ac  = subs(E,{A,X0,B,C},{Ac_1,Ac_2,Ac_3,Ac_4});

B2 = d2fz0*(EcF-EcP)/Ac - 1;
dEc = Ac*fz*(1+B2*r^4)/d2fz0;
Ec = EcP + dEc;
Fc = rho*Ec;

data.functional = subs(Fc,r,z); 
data.functional_a0 = subs(Fc,r,1);  
data.functional_b0 = subs(Fc,r,1); 
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'VWN5_C';
data.citation = 'S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980';
data.description = 'VWN5 Correlation Functional';
