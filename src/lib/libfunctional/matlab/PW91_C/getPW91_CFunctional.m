function data = getPW91_CFunctional()

data = getDefaultFunctional();
data.deriv2 = 0; %We're not ready yet

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'two_13', 'k', 'pi_m12', 'd2fz0'...
        'Aa',  'a1a', 'b1a', 'b2a', 'b3a', 'b4a',...
        'c0p', 'a1p', 'b1p', 'b2p', 'b3p', 'b4p',...
        'c0f', 'a1f', 'b1f', 'b2f', 'b3f', 'b4f',...
        'alph', 'bet', 'nu','Cc0', 'Cx', 'Cc1', 'Cc2', 'Cc3', 'Cc4', 'Cc5', 'Cc6' };  
data.param_vals = [1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), 3.0^(1.0/3.0)*pi^(2.0/3.0), ... 
    sqrt(1.0/pi), 1.709921,...
    0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, ...
    0.031091, 0.21370, 7.5957, 3.5876, 1.6382,  0.49294, ...
    0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517,...
    0.09, 16*(3*pi*pi)^(1.0/3.0)/pi*0.004235, 16*(3*pi*pi)^(1.0/3.0)/pi, 0.004235, -0.001667, 2.568, 23.266, 0.007389, 8.723, 0.472, 0.07389];
    
syms r c two_13 k pi_m12 d2fz0 real
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

Ec = EcP - Ac*fz*(1-r^4)/d2fz0 + (EcF - EcP) * fz * r^4;

% H correction
syms alph bet Cc0 nu Cx Cc1 Cc2 Cc3 Cc4 Cc5 Cc6 real
kF = k*rho^(1/3);
ks = 2*pi_m12*sqrt(kF);

gs = 1/2*((1+r)^(2/3) + (1-r)^(2/3));
T = 1/2*sqrt(gamm)/(gs*ks*rho);
A = 2*alph/bet/(exp(-2*alph*Ec/(gs^3*bet^2))-1);
Cc = 1/1000*(Cc1 + Cc2*rs+Cc3*rs^2)/(1+Cc4*rs+Cc5*rs^2+Cc6*rs^3)-Cx;
H0 = 1/2*gs^3*bet^2*log(1+2*(alph*(T^2 + A*T^4))/(bet*(1+A*T^2+A^2*T^4)))/alph;
H1 = nu*(Cc-Cc0-3/7*Cx)*gs^3*T^2*exp(-100*gs^4*ks^2*T^2/kF^2);

H = H0+H1;
pw91c = rho*(Ec+H); 

data.functional = subs(pw91c,r,z); 
data.functional_a0 = subs(pw91c,r,1);
data.functional_b0 = subs(pw91c,r,1);
data.functional_a0b0 = 0; 

data.type = 'c';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'PW91_C';
data.citation = 'J.P. Perdew, et. al., Phys. Rev. B., 46(11), 6671-6687, 1992';
data.description = 'PW91 Correlation';
