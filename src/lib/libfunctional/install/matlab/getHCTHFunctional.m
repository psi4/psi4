function data = getHCTHFunctional()

data = getDefaultFunctional();

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c0', 'two_13', 'd2fz0', 'c'...
        'Aa',  'a1a', 'b1a', 'b2a', 'b3a', 'b4a',...
        'c0p', 'a1p', 'b1p', 'b2p', 'b3p', 'b4p',...
        'c0f', 'a1f', 'b1f', 'b2f', 'b3f', 'b4f',...
        'gcab', 'gcaa', 'gx',...
        'ccab0', 'ccab1', 'ccab2', 'ccab3', 'ccab4',...
        'ccaa0', 'ccaa1', 'ccaa2', 'ccaa3', 'ccaa4',...
        'cx0', 'cx1', 'cx2', 'cx3', 'cx4'};
data.param_vals = [-3.0/8.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), 1.709921, 1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0),...
    0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, ...
    0.031091, 0.21370, 7.5957, 3.5876, 1.6382,  0.49294, ...
    0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517,...
    0.006, 0.2, 0.004,...
    0.729974, 3.35287, -11.5430, 8.08564, -4.47857,...
    0.222601, -0.0338622, -0.0125170, -0.802496, 1.55396,...
    1.09320, -0.744056, 5.59920, -6.78549, 4.49357];
syms c two_13 c0 d2fz0 real
syms zz ra rb real 

z = (rho_a - rho_b) / (rho_a + rho_b);
fz = ((1+zz)^(4/3) + (1-zz)^(4/3) -2)/(2*two_13-2);
rs = c*(ra + rb)^(-1/3);

% Ex portion
ExaLSDA = c0*rho_a^(4/3);
ExbLSDA = c0*rho_b^(4/3);

% Ec portion
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

syms cx0 cx1 cx2 cx3 cx4 gx real; 
Gxaa = subs(G,{C0, C1, C2, C3, C4, g, t}, {cx0, cx1, cx2, cx3, cx4, gx, gamma_aa/rho_a^(8/3)});
Gxbb = subs(G,{C0, C1, C2, C3, C4, g, t}, {cx0, cx1, cx2, cx3, cx4, gx, gamma_bb/rho_b^(8/3)});

% Synthesis
K = EcabLSDA*Gcab + EcaaLSDA*Gcaa + EcbbLSDA*Gcbb + ExaLSDA*Gxaa + ExbLSDA*Gxbb;
K_a0 = EcbbLSDA*Gcbb + ExbLSDA*Gxbb;  
K_b0 = EcaaLSDA*Gcaa + ExaLSDA*Gxaa;  

data.functional = subs(K, zz, z,0); 
data.functional_a0 = subs(K_a0, zz, 1, 0); 
data.functional_b0 = subs(K_b0, zz, 1, 0); 
data.functional_a0b0 = 0; 

data.type = 'xc';
data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'HCTH';
data.citation = 'F.A. Hamprecht et. al., J. Chem. Phys., 109(15), 6264-6271';
data.description = 'HCTH Power Series GGA';
