function data = getSSuperFunctional()

data = getDefaultSuperFunctional();

data.functionals = {'S'};
data.weights = [1.0];

data.type = 'x';

data.is_lsda = 1;
data.is_gga = 0;
data.is_meta = 0;
data.exact_exchange = 0.0;
data.pt2 = 0.0;
data.omega = 0.0;
data.dashD_weight = 0.0;
data.dashD = '';

data.name = 'S';
data.citation = 'J.C. Slater, Phys. Rev., 81(3):385-390, 1951';
data.description = 'Simple Slater LSDA exchange functional';
