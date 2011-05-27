function data = getM05_2XFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'X', 'C_k', 'k','e','kp','mu_', ...
    'a_1','a_2','a_3','a_4','a_5','a_6','a_7','a_8','a_9','a_10','a_11', ...
    'two_13', 'd2fz0', 'c'...
    'Aa',  'a1a', 'b1a', 'b2a', 'b3a', 'b4a',...
    'c0p', 'a1p', 'b1p', 'b2p', 'b3p', 'b4p',...
    'c0f', 'a1f', 'b1f', 'b2f', 'b3f', 'b4f',...
    'gcab', 'gcaa', ...
    'ccab0', 'ccab1', 'ccab2', 'ccab3', 'ccab4',...
    'ccaa0', 'ccaa1', 'ccaa2', 'ccaa3', 'ccaa4'};
data.param_vals = [0.56, 3.0/5.0*(6.0*pi*pi)^(2.0/3.0),3.0^(1.0/3.0)*pi^(2.0/3.0), -3.0/4.0/pi, 0.804, 0.2195149727645171,... 
    -0.56833, -1.30057, 5.50070, 9.06402, -32.21075, -23.73298, 70.22996, 29.88614, -60.25778, -13.22205, 15.23696,...
    2.0^(1.0/3.0), 1.709921, 1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0),... %Correlation
    0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, ...
    0.031091, 0.21370, 7.5957, 3.5876, 1.6382,  0.49294, ...
    0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, ...
    0.0031, 0.06, ...
    1.0, 1.09297, -3.79171, 2.82810, -10.58909,...
    1.0,-3.05430,  7.61854, 1.47665, -11.92365];
syms k e kp mu_ real;
syms n s real;

% EXCHANGE
% Start with PBE_X
kF = k*(2*n)^(1/3);
eX = e*kF;
S = 1/2*sqrt(4*s)/(kF*2*n);

PBE_X = 1 + kp - kp*(1 + mu_*S^2/kp)^-1;
PBE_a = rho_a*subs(PBE_X,{n,s},{rho_a, gamma_aa});
PBE_b = rho_b*subs(PBE_X,{n,s},{rho_b, gamma_bb});

% Now correct with Becke's unreasonably large power series
syms C_k real
tau_LDA_a = C_k*rho_a^(5/3);
tau_LDA_b = C_k*rho_b^(5/3);

t_a = tau_LDA_a/tau_a; 
t_b = tau_LDA_b/tau_b; 

w_a = (t_a - 1)/(t_a + 1);
w_b = (t_b - 1)/(t_b + 1);

syms a_1 a_2 a_3 a_4 a_5 a_6 a_7 a_8 a_9 a_10 a_11 real
H_a = 1 + a_1 * w_a^1 + a_2 * w_a^2 + a_3 * w_a^3 + a_4 * w_a^4 + a_5 * w_a^5 + a_6 * w_a^6+...
          a_7 * w_a^7 + a_8 * w_a^8 + a_9 * w_a^9 + a_10 * w_a^10 + a_11 * w_a^11;
H_b = 1 + a_1 * w_b^1 + a_2 * w_b^2 + a_3 * w_b^3 + a_4 * w_b^4 + a_5 * w_b^5 + a_6 * w_b^6+...
          a_7 * w_b^7 + a_8 * w_b^8 + a_9 * w_b^9 + a_10 * w_b^10 + a_11 * w_b^11;

M05_Xa = PBE_a*H_a;
M05_Xb = PBE_b*H_b;

% CORRELATION
syms c two_13 d2fz0 real
syms zz ra rb real 

z = (rho_a - rho_b) / (rho_a + rho_b);
fz = ((1+zz)^(4/3) + (1-zz)^(4/3) -2)/(2*two_13-2);
rs = c*(ra + rb)^(-1/3);

syms A a1 b1 b2 b3 b4 real
Gc = -2*A*(1+a1*rs)*log(1+1/(2*A*(b1*sqrt(rs)+b2*rs+b3*rs^(3/2)+b4*rs^2)));

syms Aa a1a b1a b2a b3a b4a real
Ac =  -subs(Gc, {A, a1, b1, b2, b3, b4}, {Aa, a1a, b1a, b2a, b3a, b4a});

syms c0p a1p b1p b2p b3p b4p real
EcP =  subs(Gc, {A, a1, b1, b2, b3, b4}, {c0p, a1p, b1p, b2p, b3p, b4p});

syms c0f a1f b1f b2f b3f b4f real
EcF =  subs(Gc, {A, a1, b1, b2, b3, b4}, {c0f, a1f, b1f, b2f, b3f, b4f});

EcLSDA = (ra + rb)*(EcP + Ac*fz*(1-zz^4)/d2fz0 + (EcF - EcP) * fz * zz^4);

EcabLSDA = subs(EcLSDA, {ra, rb}, {rho_a, rho_b}) - ...
           subs(subs(EcLSDA, {ra, rb}, {rho_a, 0    }), zz, 1) - ... 
           subs(subs(EcLSDA, {ra, rb}, {0    , rho_b}), zz, -1); 

EcaaLSDA = subs(subs(EcLSDA, {ra, rb}, {rho_a, 0    }), zz, 1);
EcbbLSDA = subs(subs(EcLSDA, {ra, rb}, {0    , rho_b}), zz, -1);

% Becke power series
syms C0 C1 C2 C3 C4 g t
G = C0 + C1*g*t/(1+g*t) + C2*g^2*t^2/(1+2*g*t+g^2*t^2) + C3*g^3*t^3/(1+3*g*t + 3*g^2*t^2 + g^3*t^3) + C4*g^4*t^4/ ...
    (1+4*g*t+6*g^2*t^2+4*g^3*t^3+g^4*t^4);

syms ccab0 ccab1 ccab2 ccab3 ccab4 gcab real; 
Gcab = subs(G,{C0, C1, C2, C3, C4, g, t}, {ccab0, ccab1, ccab2, ccab3, ccab4, gcab, 1/2*gamma_aa/rho_a^(8/3) + 1/2*gamma_bb/rho_b^(8/3)});

syms ccaa0 ccaa1 ccaa2 ccaa3 ccaa4 gcaa real; 
Gcaa = subs(G,{C0, C1, C2, C3, C4, g, t}, {ccaa0, ccaa1, ccaa2, ccaa3, ccaa4, gcaa, gamma_aa/rho_a^(8/3)});
Gcbb = subs(G,{C0, C1, C2, C3, C4, g, t}, {ccaa0, ccaa1, ccaa2, ccaa3, ccaa4, gcaa, gamma_bb/rho_b^(8/3)});

% D_a/2*tau_a self-interaction correction
D_a = 2*(tau_a - gamma_aa/(8*rho_a));
D_b = 2*(tau_b - gamma_bb/(8*rho_b));
F_aa = heaviside(tau_a-1E-20)*(D_a/(2*tau_a) -1) + 1;
F_bb = heaviside(tau_b-1E-20)*(D_b/(2*tau_b) -1) + 1;

% Synthesis
K = subs(EcabLSDA*Gcab + EcaaLSDA*Gcaa*F_aa + EcbbLSDA*Gcbb*F_bb,zz,z,0);
Kb = subs(EcbbLSDA*Gcbb*F_bb,zz,1,0);  
Ka = subs(EcaaLSDA*Gcaa*F_aa,zz,1,0);  

syms X
data.functional = (1-X)*(M05_Xa + M05_Xb) + K;
data.functional_a0 = (1-X)*M05_Xb + Kb;
data.functional_b0 = (1-X)*M05_Xa + Ka;
data.functional_a0b0 = 0.0; 

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 1;
data.type = 'xc';

data.name = 'M05_2X';
data.citation = 'Zhao, Y., Schultz, N. E., Truhlar, D. G., J. Chem. Theory Comput. 2, 364, 2006';
data.description = 'M05-2X Meta-GGA Functional';
