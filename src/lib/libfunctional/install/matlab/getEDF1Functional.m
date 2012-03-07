function data = getEDF1Functional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

% EDF1 params 
data.param_names{1} = 'c0';
data.param_vals(1) = -3/8*3^(1/3)*4^(2/3)*pi^(-1/3);
data.param_names{2} = 'd1';
data.param_vals(2) = .0035;
data.param_names{3} = 'd2';
data.param_vals(3) = .0042;
data.param_names{4} = 'e0';
data.param_vals(4) = 1.030952;
data.param_names{5} = 'e1';
data.param_vals(5) = 10.4017;
data.param_names{6} = 'e2';
data.param_vals(6) = -8.44793;

% EDF1 Exchange (based on B)
syms d n s real
ldax = n^(4/3);
chi = sqrt(s)/n^(4/3);
b88x = -d*chi^2/(1+6*d*chi*asinh(chi));

syms e0 c0 e1 d1 e2 d2 real 
edf1x = (e0*c0 + e1*subs(b88x,d,d1) + e2*subs(b88x,d,d2))*ldax;

% Alternate LYP Parameters
data.param_names{7} = 'A';
data.param_vals(7) = 0.055;
data.param_names{8} = 'B';
data.param_vals(8) = 0.158;
data.param_names{9} = 'C';
data.param_vals(9) = 0.25;
data.param_names{10} = 'Dd';
data.param_vals(10) = 0.3505;
data.param_names{11} = 'CFext';
data.param_vals(11) = 8.0*2.0^(2.0/3.0)*(3.0/10.0)*3.0^(2.0/3.0)*(pi*pi)^(2.0/3.0);

% LYP (with alternate parameters
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

data.functional = lyp + subs(edf1x,{n,s},{rho_a,gamma_aa})+ subs(edf1x,{n,s},{rho_b,gamma_bb}); 
data.functional_a0 = subs(edf1x,{n,s},{rho_b,gamma_bb}); 
data.functional_b0 = subs(edf1x,{n,s},{rho_a,gamma_aa});  
data.functional_a0b0 = 0; 

data.type = 'xc';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'EDF1';
data.citation = 'R.D. Adamson, P.M.W. Gill, and J.A. Pople, Chem. Phys. Lett., 284, 6-11, 1998';
data.description = 'Empirical Density Function #1';
