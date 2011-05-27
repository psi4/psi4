function data = getEDF1_CFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

% Alternate LYP Parameters
data.param_names{1} = 'A';
data.param_vals(1) = 0.055;
data.param_names{2} = 'B';
data.param_vals(2) = 0.158;
data.param_names{3} = 'C';
data.param_vals(3) = 0.25;
data.param_names{4} = 'Dd';
data.param_vals(4) = 0.3505;
data.param_names{5} = 'CFext';
data.param_vals(5) = 8.0*2.0^(2.0/3.0)*(3.0/10.0)*3.0^(2.0/3.0)*(pi*pi)^(2.0/3.0);

% LYP (with alternate parameters)
syms A B C Dd CFext omega delta rho gamm real;
lyp = ...
    -4.0*A*rho_a*rho_b/rho/(1.0+Dd/rho^(1.0/3.0)) - A*B*omega*(rho_a*rho_b*(CFext*(rho_a^(8.0/3.0)+rho_b^(8.0/3.0))...
    + (47.0/18.0 - 7.0/18.0*delta)*gamm... %GGA starts here
    -(5.0/2.0-1.0/18.0*delta)*(gamma_aa+gamma_bb) - 1.0/9.0*(delta - 11.0)*(rho_a*gamma_aa + rho_b*gamma_bb)/rho)...
    -2.0/3.0*rho*rho*gamm + (2.0/3.0*rho*rho - rho_a*rho_a)*gamma_bb + (2.0/3.0*rho*rho - rho_b*rho_b)*gamma_aa);
lyp = subs(lyp,omega,exp(-C/rho^(1.0/3.0))/(1.0+Dd/rho^(1.0/3.0))*rho^(-11.0/3.0),0);
lyp = subs(lyp,delta,C/rho^(1.0/3.0)+Dd/rho^(1.0/3.0)/(1.0+Dd/rho^(1.0/3.0)),0);
lyp = subs(lyp,rho,rho_a+rho_b,0);
lyp = subs(lyp,gamm,gamma_aa+gamma_bb+2.0*gamma_ab,0);

data.functional = lyp;
data.functional_a0 = 0; 
data.functional_b0 = 0; 
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'EDF1_C';
data.citation = 'R.D. Adamson, P.M.W. Gill, and J.A. Pople, Chem. Phys. Lett., 284, 6-11, 1998';
data.description = 'Empirical Density Function #1 (Correlation only)';
