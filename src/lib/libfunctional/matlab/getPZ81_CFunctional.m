function data = getPZ81_CFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'two_13', 'EcPld_1', 'EcPld_2','EcPld_3','EcFld_1','EcFld_2','EcFld_3',...
    'EcPhd_1', 'EcPhd_2','EcPhd_3','EcPhd_4', 'EcFhd_1','EcFhd_2','EcFhd_3','EcFhd_4'};
data.param_vals = [1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), ...
    -0.1423, 1.0529, 0.3334,...
    -0.0843, 1.3981, 0.2611,...
     0.0311, -0.048, 0.0020, -0.0116,...
     0.01555, -0.0269, 0.0007, -0.0048];

syms rho c r two_13 real
rho = rho_a + rho_b;
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

data.functional = rho*subs(Ec,r,z);
data.functional_a0 = rho*(subs(Ec,r,1)); 
data.functional_b0 = rho*(subs(Ec,r,1));  
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'PZ81_C';
data.citation = 'J.P. Perdew, A. Zunger, Phys. Rev. B., 23, 5048-5079, 1981';
data.description = 'PZ81 Correlation';
