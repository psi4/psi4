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

% Symbolic tree derivativesi
expand = data.expand;
root = 'functional';
params = data.param_names;

%========== ORDER 0 ==========%

array = 'functional_[index]';
buffer{1} = codeTree(data.functional, root, array, expand, true);
buffer{2} = codeTree(data.functional_a0, root, array, expand);
buffer{3} = codeTree(data.functional_b0, root, array, expand);
replaceInFile('UKS_FUNCTIONAL',buildUKSII(buffer),[name '_functional.cc']);

%========== ORDER 1 ==========%

partial1 = 'rho_a';
array = 'v_rho_a_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_RHO_A',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
array = 'v_rho_b_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_RHO_B',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
array = 'v_gamma_aa_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_GAMMA_AA',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_ab';
array = 'v_gamma_ab_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_GAMMA_AB',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_bb';
array = 'v_gamma_bb_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_GAMMA_BB',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'tau_a';
array = 'v_tau_a_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_TAU_A',buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'tau_b';
array = 'v_tau_b_[index]';
v = order1Tree(data.functional, root, partial1, params);
v_a0 = order1Tree(data.functional_a0, root, partial1, params);
v_b0 = order1Tree(data.functional_b0, root, partial1, params);
buffer{1} = codeTree(v, [root '_' partial1], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1], array, expand);
replaceInFile('UKS_V1_TAU_B',buildUKSII(buffer),[name '_functional.cc']);

%========== ORDER 2 ==========%

partial1 = 'rho_a';
partial2 = 'rho_a';
array = 'v_rho_a_rho_a_[index]';
marker = 'UKS_V2_RHO_A_RHO_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'rho_b';
array = 'v_rho_a_rho_b_[index]';
marker = 'UKS_V2_RHO_A_RHO_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'rho_b';
array = 'v_rho_b_rho_b_[index]';
marker = 'UKS_V2_RHO_B_RHO_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'gamma_aa';
array = 'v_rho_a_gamma_aa_[index]';
marker = 'UKS_V2_RHO_A_GAMMA_AA';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'gamma_ab';
array = 'v_rho_a_gamma_ab_[index]';
marker = 'UKS_V2_RHO_A_GAMMA_AB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'gamma_bb';
array = 'v_rho_a_gamma_bb_[index]';
marker = 'UKS_V2_RHO_A_GAMMA_BB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'gamma_aa';
array = 'v_rho_b_gamma_aa_[index]';
marker = 'UKS_V2_RHO_B_GAMMA_AA';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'gamma_ab';
array = 'v_rho_b_gamma_ab_[index]';
marker = 'UKS_V2_RHO_B_GAMMA_AB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'gamma_bb';
array = 'v_rho_b_gamma_bb_[index]';
marker = 'UKS_V2_RHO_B_GAMMA_BB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
partial2 = 'gamma_aa';
array = 'v_gamma_aa_gamma_aa_[index]';
marker = 'UKS_V2_GAMMA_AA_GAMMA_AA';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
partial2 = 'gamma_ab';
array = 'v_gamma_aa_gamma_ab_[index]';
marker = 'UKS_V2_GAMMA_AA_GAMMA_AB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
partial2 = 'gamma_bb';
array = 'v_gamma_aa_gamma_bb_[index]';
marker = 'UKS_V2_GAMMA_AA_GAMMA_BB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_ab';
partial2 = 'gamma_ab';
array = 'v_gamma_ab_gamma_ab_[index]';
marker = 'UKS_V2_GAMMA_AB_GAMMA_AB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_ab';
partial2 = 'gamma_bb';
array = 'v_gamma_ab_gamma_bb_[index]';
marker = 'UKS_V2_GAMMA_AB_GAMMA_BB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_bb';
partial2 = 'gamma_bb';
array = 'v_gamma_bb_gamma_bb_[index]';
marker = 'UKS_V2_GAMMA_BB_GAMMA_BB';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'tau_a';
array = 'v_rho_a_tau_a_[index]';
marker = 'UKS_V2_RHO_A_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_a';
partial2 = 'tau_b';
array = 'v_rho_a_tau_b_[index]';
marker = 'UKS_V2_RHO_A_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'tau_a';
array = 'v_rho_b_tau_a_[index]';
marker = 'UKS_V2_RHO_B_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'rho_b';
partial2 = 'tau_b';
array = 'v_rho_b_tau_b_[index]';
marker = 'UKS_V2_RHO_B_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'tau_a';
partial2 = 'tau_a';
array = 'v_tau_a_tau_a_[index]';
marker = 'UKS_V2_TAU_A_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'tau_a';
partial2 = 'tau_b';
array = 'v_tau_a_tau_b_[index]';
marker = 'UKS_V2_TAU_A_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'tau_b';
partial2 = 'tau_b';
array = 'v_tau_b_tau_b_[index]';
marker = 'UKS_V2_TAU_B_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
partial2 = 'tau_a';
array = 'v_gamma_aa_tau_a_[index]';
marker = 'UKS_V2_GAMMA_AA_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_ab';
partial2 = 'tau_a';
array = 'v_gamma_ab_tau_a_[index]';
marker = 'UKS_V2_GAMMA_AB_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_bb';
partial2 = 'tau_a';
array = 'v_gamma_bb_tau_a_[index]';
marker = 'UKS_V2_GAMMA_BB_TAU_A';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_aa';
partial2 = 'tau_b';
array = 'v_gamma_aa_tau_b_[index]';
marker = 'UKS_V2_GAMMA_AA_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_ab';
partial2 = 'tau_b';
array = 'v_gamma_ab_tau_a_[index]';
marker = 'UKS_V2_GAMMA_AB_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

partial1 = 'gamma_bb';
partial2 = 'tau_b';
array = 'v_gamma_bb_tau_b_[index]';
marker = 'UKS_V2_GAMMA_BB_TAU_B';
v = order2Tree(data.functional, root, partial1, partial2, params);
v_a0 = order2Tree(data.functional_a0, root, partial1, partial2, params);
v_b0 = order2Tree(data.functional_b0, root, partial1, partial2, params);
buffer{1} = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
buffer{2} = codeTree(v_a0, [root '_' partial1 '_' partial2], array, expand);
buffer{3} = codeTree(v_b0, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildUKSII(buffer),[name '_functional.cc']);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       RKS FUNCTIONALS
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

expand = data.expand;
root = 'functional';
params = data.param_names;

%========== ORDER 0 ==========%

% zk <- zk
array = 'functional_[index]';
marker = 'RKS_FUNCTIONAL';
v = cleanRKSTree(data.functional, [root], params);
buffer = codeTree(v, [root], array, expand, true);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

%========== ORDER 1 ==========%

% v_rho_a <- v_rho_a
partial1 = 'rho_a';
array = 'v_rho_a_[index]';
marker = 'RKS_V1_RHO_A';
v = order1Tree(data.functional, root, partial1, params);
v = cleanRKSTree(v, [root '_' partial1], params);
buffer = codeTree(v, [root '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_gamma_aa <- 2.0 * v_gamma_aa + v_gamma_ab
partial1 = 'gamma_aa';
partial2 = 'gamma_ab';
array = 'v_gamma_aa_[index]';
marker = 'RKS_V1_GAMMA_A';
v1 = order1Tree(data.functional, root, partial1, params);
v2 = order1Tree(data.functional, root, partial2, params);
v = addTrees({v1 v2},{'functional_gamma_aa', 'functional_gamma_ab'}, ...
        [root '_' partial1], [2 1]); 
v = cleanRKSTree(v, [root '_' partial1], params);
buffer = codeTree(v, [root '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_tau_a <- v_tau_a
partial1 = 'tau_a';
array = 'v_tau_a_[index]';
marker = 'RKS_V1_TAU_A';
v = order1Tree(data.functional, root, partial1, params);
v = cleanRKSTree(v, [root '_' partial1], params);
buffer = codeTree(v, [root '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

%========== ORDER 2 ==========%

% v_rho_a_rho_a <- v_rho_a_rho_a + v_rho_a_rho_b 
partial1 = 'rho_a';
partial2 = 'rho_b';
array = 'v_rho_a_rho_a_[index]';
marker = 'RKS_V2_RHO_A_RHO_A';
v1 = order2Tree(data.functional, root, partial1, partial1, params);
v2 = order2Tree(data.functional, root, partial1, partial2, params);
v = addTrees({v1 v2},{'functional_rho_a_rho_a', 'functional_rho_a_rho_b'}, ...
        [root '_' partial1 '_' partial1], [1 1]);
v = cleanRKSTree(v, [root '_' partial1 '_' partial1], params);
buffer = codeTree(v, [root '_' partial1 '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_rho_a_gamma_aa <-  v_rho_a_gamma_aa + v_rho_a_gamma_ab + v_rho_a_gamma_bb 
partial1 = 'rho_a';
partial2 = 'gamma_aa';
partial3 = 'gamma_ab';
partial4 = 'gamma_bb';
array = 'v_rho_a_gamma_aa_[index]';
marker = 'RKS_V2_RHO_A_GAMMA_AA';
v1 = order2Tree(data.functional, root, partial1, partial2, params);
v2 = order2Tree(data.functional, root, partial1, partial3, params);
v3 = order2Tree(data.functional, root, partial1, partial4, params);
v = addTrees({v1 v2 v3},{'functional_rho_a_gamma_aa', 'functional_rho_a_gamma_ab', 'functional_rho_a_gamma_bb'}, ...
        [root '_' partial1 '_' partial2], [1 1 1]);
v = cleanRKSTree(v, [root '_' partial1 '_' partial2], params);
buffer = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_gamma_aa_gamma_aa <- 2 v_gamma_aa_gamma_aa + 4 v_gamma_aa_gamma_ab + 2 v_gamma_aa_gamma_bb + v_gamma_ab_gamma_ab
partial1 = 'gamma_aa';
partial2 = 'gamma_ab';
partial3 = 'gamma_bb';
array = 'v_gamma_aa_gamma_aa_[index]';
marker = 'RKS_V2_GAMMA_AA_GAMMA_AA';
v1 = order2Tree(data.functional, root, partial1, partial1, params);
v2 = order2Tree(data.functional, root, partial1, partial2, params);
v3 = order2Tree(data.functional, root, partial1, partial3, params);
v4 = order2Tree(data.functional, root, partial2, partial2, params);
v = addTrees({v1 v2 v3 v4},{'functional_gamma_aa_gamma_aa', 'functional_gamma_aa_gamma_ab',...
         'functional_gamma_aa_gamma_bb', 'functional_gamma_ab_gamma_ab'}, ...
        [root '_' partial1 '_' partial1], [2 4 2 1]);
v = cleanRKSTree(v, [root '_' partial1 '_' partial1], params);
buffer = codeTree(v, [root '_' partial1 '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_rho_a_tau_a <- v_rho_a_tau_a + v_rho_a_tau_b (I think)
partial1 = 'rho_a';
partial2 = 'tau_a';
partial3 = 'tau_b';
array = 'v_rho_a_tau_a_[index]';
marker = 'RKS_V2_RHO_A_TAU_A';
v1 = order2Tree(data.functional, root, partial1, partial2, params);
v2 = order2Tree(data.functional, root, partial1, partial3, params);
v = addTrees({v1 v2},{'functional_rho_a_tau_a', 'functional_rho_a_tau_b'}, ...
        [root '_' partial1 '_' partial2], [1 1]);
v = cleanRKSTree(v, [root '_' partial1 '_' partial2], params);
buffer = codeTree(v, [root '_' partial1 '_' partial2], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_tau_a_tau_a <- v_tau_a_tau_a + v_tau_a_tau_b 
partial1 = 'tau_a';
partial2 = 'tau_b';
array = 'v_tau_a_tau_a_[index]';
marker = 'RKS_V2_TAU_A_TAU_A';
v1 = order2Tree(data.functional, root, partial1, partial1, params);
v2 = order2Tree(data.functional, root, partial1, partial2, params);
v = addTrees({v1 v2},{'functional_tau_a_tau_a', 'functional_tau_a_tau_b'}, ...
        [root '_' partial1 '_' partial1], [1 1]);
v = cleanRKSTree(v, [root '_' partial1 '_' partial1], params);
buffer = codeTree(v, [root '_' partial1 '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);

% v_tau_a_gamma_aa <-  v_tau_a_gamma_aa + v_tau_a_gamma_ab + v_tau_a_gamma_bb 
partial1 = 'tau_a';
partial2 = 'gamma_aa';
partial3 = 'gamma_ab';
partial4 = 'gamma_bb';
array = 'v_gamma_aa_tau_a_[index]';
marker = 'RKS_V2_GAMMA_AA_TAU_A';
v1 = order2Tree(data.functional, root, partial1, partial2, params);
v2 = order2Tree(data.functional, root, partial1, partial3, params);
v3 = order2Tree(data.functional, root, partial1, partial4, params);
v = addTrees({v1 v2 v3},{'functional_tau_a_gamma_aa', 'functional_tau_a_gamma_ab', 'functional_tau_a_gamma_bb'}, ...
        [root '_' partial2 '_' partial1], [1 1 1]);
v = cleanRKSTree(v, [root '_' partial2 '_' partial1], params);
buffer = codeTree(v, [root '_' partial2 '_' partial1], array, expand);
replaceInFile(marker, buildRKSII(buffer),[name '_functional.cc']);


