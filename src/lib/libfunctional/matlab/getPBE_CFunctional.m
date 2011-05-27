function data = getPBE_CFunctional()

data = getDefaultFunctional();
data.deriv2 = 0; %We're not ready yet

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'two_13', 'k', 'pi_m12', 'bet', 'gammas', 'd2fz0'...
        'Aa', 'a1a', 'b1a', 'b2a', 'b3a', 'b4a',...
        'c0p', 'a1p', 'b1p', 'b2p', 'b3p', 'b4p',...
        'c0f', 'a1f', 'b1f', 'b2f', 'b3f', 'b4f'};  
data.param_vals = [1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), 3.0^(1.0/3.0)*pi^(2.0/3.0), ... 
    sqrt(1.0/pi), 0.0667245506031492, (1.0-log(2.0))/(pi*pi),4.0/9.0/(2.0^(1.0/3.0)-1.0),...
    0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, ...
    0.0310907, 0.21370, 7.5957, 3.5876, 1.6382,  0.49294, ...
    0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517];
    
syms r c two_13 k pi_m12 bet gammas d2fz0 real
rho = rho_a + rho_b;
gamm = gamma_aa + gamma_bb + 2*gamma_ab;

z = (rho_a-rho_b)/rho;
fz = ((1+r)^(4/3) + (1-r)^(4/3) -2)/(2*two_13-2);
rs = c*rho^(-1/3);

% Ec portion
syms A a1 b1 b2 b3 b4 real
Gc = -2*A*(1+a1*rs)*log(1+1/(2*A*(b1*sqrt(rs)+b2*rs+b3*rs^(3/2)+b4*rs^2)));

syms Aa a1a b1a b2a b3a b4a real
Ac =  subs(Gc, {A, a1, b1, b2, b3, b4}, {Aa, a1a, b1a, b2a, b3a, b4a});

syms c0p a1p b1p b2p b3p b4p real
EcP =  subs(Gc, {A, a1, b1, b2, b3, b4}, {c0p, a1p, b1p, b2p, b3p, b4p});

syms c0f a1f b1f b2f b3f b4f real
EcF =  subs(Gc, {A, a1, b1, b2, b3, b4}, {c0f, a1f, b1f, b2f, b3f, b4f});

Ec = EcP + Ac*fz*(1-r^4)/d2fz0 + (EcF - EcP) * fz * r^4;

% H correction
kF = k*rho^(1/3);
ks = 2*pi_m12*sqrt(kF);

phi = 1/2*((1+r)^(2/3) + (1-r)^(2/3));
T = 1/2*sqrt(gamm)/(phi*ks*rho);
A = bet/gammas/(exp(-Ec/(gammas*phi^3)) - 1);
H = gammas*phi^3*log(1 + (bet*T^2*(1 + A*T^2))/(gammas*(1 + A*T^2 + A^2*T^4)));
Q = Ec + H;

data.functional = rho*subs(Q,r,z); 
data.functional_a0 = rho*subs(Q,{r, rho_a, gamma_aa},{1, 0, 0});
data.functional_b0 = rho*subs(Q,{r, rho_b, gamma_bb},{1, 0, 0});
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'PBE_C';
data.citation = 'J.P. Perdew, K. Burke., and M. Ernzerhof, Phys. Rev. Lett., 77, 3865-3868, 1996';
data.description = 'PBE Correlation';
