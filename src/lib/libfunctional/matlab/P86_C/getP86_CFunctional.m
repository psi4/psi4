function data = getP86_CFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'two_13', 'EcPld_1', 'EcPld_2','EcPld_3','EcFld_1','EcFld_2','EcFld_3',...
    'EcPhd_1', 'EcPhd_2','EcPhd_3','EcPhd_4', 'EcFhd_1','EcFhd_2','EcFhd_3','EcFhd_4', ...
    'Fg', 'Bg', 'Cx', 'Cinf', 'Cg_1', 'Cg_2', 'Cg_3','Cg_4', 'Pg_1'};
data.param_vals = [1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), ...
    -0.1423, 1.0529, 0.3334,...
    -0.0843, 1.3981, 0.2611,...
     0.0311, -0.048, 0.0020, -0.0116,...
     0.01555, -0.0269, 0.0007, -0.0048, ...
     0.11, 0.000007389, 0.001667, 0.004235,...
     0.002568, 0.023266, 8.723, 0.472, 1.745];

syms rho c r two_13 real
rho = rho_a + rho_b;
gamm = gamma_aa + gamma_bb + 2*gamma_ab;
% PZ81 LSDA
z = (rho_a-rho_b)/rho;
fz = ((1+r)^(4/3) + (1-r)^(4/3) -2)/(2*two_13-2);
rs = c*rho^(-1/3);

syms Bld Cld Dld Ahd Bhd Chd Dhd x real
Eld = Cld/(1+Bld*sqrt(rs)+Dld*rs);
Ehd = Ahd*log(rs) + Bhd + Chd*rs*log(rs) + Dhd*rs;

syms EcPld_1 EcPld_2 EcPld_3 EcFld_1 EcFld_2 EcFld_3 real
syms EcPhd_1 EcPhd_2 EcPhd_3 EcPhd_4 EcFhd_1 EcFhd_2 EcFhd_3 EcFhd_4 real
EcPld = subs(Eld,{Cld,Bld,Dld},{EcPld_1, EcPld_2, EcPld_3});
EcFld = subs(Eld,{Cld,Bld,Dld},{EcFld_1, EcFld_2, EcFld_3});
EcPhd = subs(Ehd,{Ahd,Bhd,Chd,Dhd}, {EcPhd_1, EcPhd_2, EcPhd_3, EcPhd_4});
EcFhd = subs(Ehd,{Ahd,Bhd,Chd,Dhd}, {EcFhd_1, EcFhd_2, EcFhd_3, EcFhd_4});

Ec = heaviside(1-rs)*(EcPhd + (EcFhd - EcPhd)*fz) + ...
     heaviside(rs-1)*(EcPld + (EcFld - EcPld)*fz);

% P86 GGA Correction
dz = two_13*sqrt((1/2+1/2*r)^(5/3) + (1/2-1/2*r)^(5/3));

syms Fg Bg Cx Cinf Cg_1 Cg_2 Cg_3 Cg_4 Pg_1 real
Cg = Cx + (Cg_1 + Cg_2*rs+Bg*rs^2)/(1+Cg_3*rs+Cg_4*rs^2+10000*Bg*rs^3);
Pg = Pg_1*Fg*Cinf*sqrt(gamm)/(Cg*rho^(7/6));
dP86c = exp(-Pg)*Cg*gamm/(rho^(4/3)*dz);

P86c = rho*Ec + dP86c;

data.functional = subs(P86c,r,z);
data.functional_a0 = (subs(P86c,r,1)); 
data.functional_b0 = (subs(P86c,r,1));  
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'P86_C';
data.citation = 'J.P. Perdew, Phys. Rev. B., 33, 8822-8824, 1986';
data.description = 'P86 Correlation (PZ81 LSDA + P86 GGA)';
