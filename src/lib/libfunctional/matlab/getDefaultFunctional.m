function data = getDefaultFunctional()

data.algorithm = 1;
data.expand = 0;

%Functional symbolic definition
data.functional = 0;
%Functional symbolic definition (for rho_a -> 0)
data.functional_a0 = 0;
%Functional symbolic definition (for rho_b -> 0)
data.functional_b0 = 0;
%Functional symbolic definition (for rho_b  + rho_b -> 0)
data.functional_a0b0 = 0;

%Functional symbolic variables
rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

data.param_names = {};
data.param_vals = [];

%is lsda?
data.is_lsda = 0;

%is gga?
data.is_gga = 0;

%is meta?
data.is_meta = 0;

%exchange or correlation or both?
%['x', 'c', 'xc']
data.type = 'x';

% Do second partials?
data.deriv2 = 1;

%Functional name
data.name = 'null';

%Functional citation
data.citation = 'Parrish, Nature 2010';

%Functional description
data.description = 'Null functional';
