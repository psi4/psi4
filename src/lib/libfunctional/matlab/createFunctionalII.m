function createFunctionalII(data)

% Assume we're in the right directory for creation, that template.cc
% and template.h exist

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%         MATLAB VARIABLE SETUP
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Functional symbolic variables
rho_a = sym('rho_a');
rho_b = sym('rho_b');
gamma_aa = sym('gamma_aa');
gamma_bb = sym('gamma_bb');
gamma_ab = sym('gamma_ab');
tau_a = sym('tau_a');
tau_b = sym('tau_b');

is_gga = data.is_gga;
is_meta = data.is_meta;

% Functional name
name = data.name;

% Functional citation
citation = data.citation;

% Functional citation
description = data.description;

%>>>>>>>>>>>>>>>>>>>>>>>
%        H File
%<<<<<<<<<<<<<<<<<<<<<<<

% Copy template to desired .cc file
current_dir = dir('./');

if (any(strcmpi([name '_functional.h'],{current_dir.name})))
    rm_str = sprintf('rm "%s_functional.h"',name);
    system(rm_str);
end

mv_str = sprintf('cp functional.h.template "%s_functional.h"',name);
system(mv_str);

% Set name, date
name_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.h"','NAME',name,name);
system(name_str);

date_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.h"','DATE',date,name);
system(date_str);

%>>>>>>>>>>>>>>>>>>>>>>>
%       CC File
%<<<<<<<<<<<<<<<<<<<<<<<

% Copy template to desired .cc file
current_dir = dir('./');

if (any(strcmpi([name '_functional.cc'],{current_dir.name})))
    rm_str = sprintf('rm "%s_functional.cc"',name);
    system(rm_str);
end

mv_str = sprintf('cp functional.cc.template "%s_functional.cc"',name);
system(mv_str);

% Set name, date, citation, description
name_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','NAME',name,name);
system(name_str);

date_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','DATE',date,name);
system(date_str);

citation_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','CITATION',citation,name);
system(citation_str);

description_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','DESCRIPTION',description,name);
system(description_str);

% Set GGA/META
if (is_gga)
    gga_val = 'true';
else
    gga_val = 'false';
end

gga_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','IS_GGA',gga_val,name);
system(gga_str);

if (is_meta)
    meta_val = 'true';
else
    meta_val = 'false';
end

meta_str = sprintf('sed -i "s/%s/%s/g" "%s_functional.cc"','IS_META',meta_val,name);
system(meta_str);

% Set definition of default parameters
params = '';
for k = 1:length(data.param_vals)
    params = [params sprintf('double %s = %22.16E;\n',data.param_names{k},data.param_vals(k))];
    params = [params sprintf('params_.push_back(make_pair("%s",%s));\n',data.param_names{k},...
        data.param_names{k})];
end

replaceInFile('DEFINE_PARAMS',params,[name '_functional.cc']);

% Set extraction of default parameters
params = '';
for k = 1:length(data.param_vals)
    params = [params sprintf('double %s = params_[%d].second;\n',data.param_names{k},k-1)];
end

replaceInFile('EXTRACT_PARAMS',params,[name '_functional.cc']);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       UKS FUNCTIONALS
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

tree = AD(data.functional, 'functional', data.param_names);
tree_a0 = AD(data.functional_a0, 'functional', data.param_names);
tree_b0 = AD(data.functional_b0, 'functional', data.param_names);

code = UKSCodeRev(tree, tree_a0, tree_b0);
replaceInFile('UKS_CODE',code,[name '_functional.cc']);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       RKS FUNCTIONALS
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%TODO

