function data = getPW92II_CFunctional()

data = getDefaultFunctional();
data.algorithm = 2;

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {'c', 'two_13', 'd2fz0'...
        'Aa', 'a1a', 'b1a', 'b2a', 'b3a', 'b4a',...
        'c0p', 'a1p', 'b1p', 'b2p', 'b3p', 'b4p',...
        'c0f', 'a1f', 'b1f', 'b2f', 'b3f', 'b4f'};  
data.param_vals = [1.0/4.0*3.0^(1.0/3.0)*4.0^(2.0/3.0)*pi^(-1.0/3.0), 2.0^(1.0/3.0), ... 
    1.709921,...    
    0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, ...
    0.031091, 0.21370, 7.5957, 3.5876, 1.6382,  0.49294, ...
    0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517];
    
syms z fz rs c two_13 d2fz0 rho  

index = 1;
tree(index).name = 'rho';
tree(index).expr = rho_a + rho_b;
index = index + 1;
tree(index).name = 'z';
tree(index).expr = (rho_a - rho_b)/(rho_a + rho_b); 
index = index + 1;
tree(index).name = 'fz';
tree(index).expr = ((1+z)^(4/3) + (1-z)^(4/3) -2)/(2*two_13-2);
index = index + 1;
tree(index).name = 'rs';
tree(index).expr = c*(rho_a + rho_b)^(-1/3); 

% Ec portion
syms rs Ac EcP EcF Ec
syms A a1 b1 b2 b3 b4 real
Gc = -2*A*(1+a1*rs)*log(1+1/(2*A*(b1*sqrt(rs)+b2*rs+b3*rs^(3/2)+b4*rs^2)));

syms Aa a1a b1a b2a b3a b4a real
syms c0p a1p b1p b2p b3p b4p real
syms c0f a1f b1f b2f b3f b4f real

index = index + 1;
tree(index).name = 'Ac';
tree(index).expr = subs(Gc, {A, a1, b1, b2, b3, b4}, {Aa, a1a, b1a, b2a, b3a, b4a}); 
index = index + 1;
tree(index).name = 'EcP';
tree(index).expr = subs(Gc, {A, a1, b1, b2, b3, b4}, {c0p, a1p, b1p, b2p, b3p, b4p}); 
index = index + 1;
tree(index).name = 'EcF';
tree(index).expr = subs(Gc, {A, a1, b1, b2, b3, b4}, {c0f, a1f, b1f, b2f, b3f, b4f});
index = index + 1;
tree(index).name = 'Ec';
tree(index).expr = EcP - Ac*fz*(1-z^4)/d2fz0 + (EcF - EcP) * fz * z^4;
index = index + 1;
tree(index).name = 'functional';
tree(index).expr = rho*Ec;

tree = buildDependencies(tree, data.param_names);

tree_a0 = tree;
tree_a0 = subsTreeElement(tree_a0, sym('z'), -1);
tree_a0 = subsTreeVariable(tree_a0,sym('rho_a'), 0);

tree_b0 = tree;
tree_b0 = subsTreeElement(tree_b0, sym('z'), 1);
tree_b0 = subsTreeVariable(tree_b0,sym('rho_b'), 0);

% Build dependencies in trees
tree = buildDependencies(tree, data.param_names);
tree_a0 = buildDependencies(tree_a0, data.param_names);
tree_b0 = buildDependencies(tree_b0, data.param_names);

% Remove orphan nodes
tree = removeOrphans(tree, 'functional');
tree_a0 = removeOrphans(tree_a0, 'functional');
tree_b0 = removeOrphans(tree_b0, 'functional');

data.functional = tree;
data.functional_a0 = tree_a0;
data.functional_b0 = tree_b0;
    
data.type = 'c';
data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.is_exchange = 0;

data.name = 'PW92II_C';
data.citation = 'J.P. Perdew and Y. Wang, Phys. Rev. B., 45(23), 13244, 1992';
data.description = 'PW92 LSDA Correlation (AD)';

