function data = getDefaultFunctional()

data.functionals = {};
data.weights = [];

data.type = 'x';

data.is_lsda = 0; 
data.is_gga = 0; 
data.is_meta = 0; 
data.exact_exchange = 0.0;
data.pt2 = 0.0;
data.omega = 0.0;
data.dashD_weight = 0.0;
data.dashD = '';

%SuperFunctional name
data.name = 'null';

%SuperFunctional citation
data.citation = 'Parrish, Nature 2010';

%SuperFunctional description
data.description = 'Null superfunctional';
