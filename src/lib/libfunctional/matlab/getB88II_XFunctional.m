function data = getB88II_XFunctional()

data = getDefaultFunctional();
data.algorithm = 2;

rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');

data.param_names{1} = 'd';
data.param_vals(1) = 0.0042;

d = sym('d');
r_a_43 = sym('r_a_43');
r_b_43 = sym('r_b_43');
chi_a = sym('chi_a');
chi_b = sym('chi_b');
asinh_a = sym('asinh_a');
asinh_b = sym('asinh_b');

index = 1;
tree(index).name = 'r_a_43';
tree(index).expr = rho_a^(4/3);
index = index + 1;
tree(index).name = 'r_b_43';
tree(index).expr = rho_b^(4/3);
index = index + 1;
tree(index).name = 'chi_a';
tree(index).expr = sqrt(gamma_aa)/rho_a^(4/3);
index = index + 1;
tree(index).name = 'chi_b';
tree(index).expr = sqrt(gamma_bb)/rho_b^(4/3);
index = index + 1;
tree(index).name = 'asinh_a';
tree(index).expr = asinh(chi_a);
index = index + 1;
tree(index).name = 'asinh_b';
tree(index).expr = asinh(chi_b);
index = index + 1;
tree(index).name = 'functional';
tree(index).expr =  -d*r_a_43*chi_a^2/(1.0+6.0*d*chi_a*asinh_a)...
                    -d*r_b_43*chi_b^2/(1.0+6.0*d*chi_b*asinh_b);

tree_a0 = tree;
tree_a0 = subsTreeElement(tree_a0, sym('functional'), -d*r_b_43*chi_b^2/(1.0+6.0*d*chi_b*asinh_b)); 

tree_b0 = tree;
tree_b0 = subsTreeElement(tree_b0, sym('functional'), -d*r_a_43*chi_a^2/(1.0+6.0*d*chi_a*asinh_a)); 

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

data.is_lsda = 1;
data.is_gga = 1;
data.is_meta = 0;
data.type = 'x';

data.name = 'B88II_X';
data.citation = 'A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988';
data.description = 'Becke Exchange GGA Correction (B88)';
