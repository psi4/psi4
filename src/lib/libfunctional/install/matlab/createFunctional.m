function createFunctional(data)

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

mv_str = sprintf('mv functional.h.template "%s_functional.h"',name);
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

mv_str = sprintf('mv functional.cc.template "%s_functional.cc"',name);
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

%========== ORDER 0 ==========%

ccode(data.functional,'file','functional');

%========== ORDER 1 ==========%

ccode(diff(data.functional,rho_a,1),'file','v_rho_a');
ccode(diff(data.functional,rho_b,1),'file','v_rho_b');
ccode(diff(data.functional,gamma_aa,1),'file','v_gamma_aa');
ccode(diff(data.functional,gamma_ab,1),'file','v_gamma_ab');
ccode(diff(data.functional,gamma_bb,1),'file','v_gamma_bb');
ccode(diff(data.functional,tau_a,1),'file','v_tau_a');
ccode(diff(data.functional,tau_b,1),'file','v_tau_b');

%========== ORDER 2 ==========%

if (data.deriv2)

    ccode(diff(diff(data.functional,rho_a,1),rho_a,1),'file','v_rho_a_rho_a');
    ccode(diff(diff(data.functional,rho_a,1),rho_b,1),'file','v_rho_a_rho_b');
    ccode(diff(diff(data.functional,rho_b,1),rho_b,1),'file','v_rho_b_rho_b');
    ccode(diff(diff(data.functional,rho_a,1),gamma_aa,1),'file','v_rho_a_gamma_aa');
    ccode(diff(diff(data.functional,rho_a,1),gamma_ab,1),'file','v_rho_a_gamma_ab');
    ccode(diff(diff(data.functional,rho_a,1),gamma_bb,1),'file','v_rho_a_gamma_bb');
    ccode(diff(diff(data.functional,rho_b,1),gamma_aa,1),'file','v_rho_b_gamma_aa');
    ccode(diff(diff(data.functional,rho_b,1),gamma_ab,1),'file','v_rho_b_gamma_ab');
    ccode(diff(diff(data.functional,rho_b,1),gamma_bb,1),'file','v_rho_b_gamma_bb');
    ccode(diff(diff(data.functional,gamma_aa,1),gamma_aa,1),'file','v_gamma_aa_gamma_aa');
    ccode(diff(diff(data.functional,gamma_aa,1),gamma_ab,1),'file','v_gamma_aa_gamma_ab');
    ccode(diff(diff(data.functional,gamma_aa,1),gamma_bb,1),'file','v_gamma_aa_gamma_bb');
    ccode(diff(diff(data.functional,gamma_ab,1),gamma_ab,1),'file','v_gamma_ab_gamma_ab');
    ccode(diff(diff(data.functional,gamma_ab,1),gamma_bb,1),'file','v_gamma_ab_gamma_bb');
    ccode(diff(diff(data.functional,gamma_bb,1),gamma_bb,1),'file','v_gamma_bb_gamma_bb');
    ccode(diff(diff(data.functional,rho_a,1),tau_a,1),'file','v_rho_a_tau_a');
    ccode(diff(diff(data.functional,rho_a,1),tau_b,1),'file','v_rho_a_tau_b');
    ccode(diff(diff(data.functional,rho_b,1),tau_a,1),'file','v_rho_b_tau_a');
    ccode(diff(diff(data.functional,rho_b,1),tau_b,1),'file','v_rho_b_tau_b');
    ccode(diff(diff(data.functional,tau_a,1),tau_a,1),'file','v_tau_a_tau_a');
    ccode(diff(diff(data.functional,tau_a,1),tau_b,1),'file','v_tau_a_tau_b');
    ccode(diff(diff(data.functional,tau_b,1),tau_b,1),'file','v_tau_b_tau_b');
    ccode(diff(diff(data.functional,gamma_aa,1),tau_a,1),'file','v_gamma_aa_tau_a');
    ccode(diff(diff(data.functional,gamma_ab,1),tau_a,1),'file','v_gamma_ab_tau_a');
    ccode(diff(diff(data.functional,gamma_bb,1),tau_a,1),'file','v_gamma_bb_tau_a');
    ccode(diff(diff(data.functional,gamma_aa,1),tau_b,1),'file','v_gamma_aa_tau_b');
    ccode(diff(diff(data.functional,gamma_ab,1),tau_b,1),'file','v_gamma_ab_tau_b');
    ccode(diff(diff(data.functional,gamma_bb,1),tau_b,1),'file','v_gamma_bb_tau_b');

end

%========== ORDER 0 ==========%

try
    ccode(data.functional_a0,'file','functional_a0');
catch exception
    syms qqqq
    ccode(diff(data.functional_a0,qqqq),'file','functional_a0')
end

%========== ORDER 1 ==========%

ccode(diff(data.functional_a0,rho_a,1),'file','v_rho_a_a0');
ccode(diff(data.functional_a0,rho_b,1),'file','v_rho_b_a0');
ccode(diff(data.functional_a0,gamma_aa,1),'file','v_gamma_aa_a0');
ccode(diff(data.functional_a0,gamma_ab,1),'file','v_gamma_ab_a0');
ccode(diff(data.functional_a0,gamma_bb,1),'file','v_gamma_bb_a0');
ccode(diff(data.functional_a0,tau_a,1),'file','v_tau_a_a0');
ccode(diff(data.functional_a0,tau_b,1),'file','v_tau_b_a0');

%========== ORDER 2 ==========%

if (data.deriv2)

    ccode(diff(diff(data.functional_a0,rho_a,1),rho_a,1),'file','v_rho_a_rho_a_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),rho_b,1),'file','v_rho_a_rho_b_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),rho_b,1),'file','v_rho_b_rho_b_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),gamma_aa,1),'file','v_rho_a_gamma_aa_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),gamma_ab,1),'file','v_rho_a_gamma_ab_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),gamma_bb,1),'file','v_rho_a_gamma_bb_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),gamma_aa,1),'file','v_rho_b_gamma_aa_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),gamma_ab,1),'file','v_rho_b_gamma_ab_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),gamma_bb,1),'file','v_rho_b_gamma_bb_a0');
    ccode(diff(diff(data.functional_a0,gamma_aa,1),gamma_aa,1),'file','v_gamma_aa_gamma_aa_a0');
    ccode(diff(diff(data.functional_a0,gamma_aa,1),gamma_ab,1),'file','v_gamma_aa_gamma_ab_a0');
    ccode(diff(diff(data.functional_a0,gamma_aa,1),gamma_bb,1),'file','v_gamma_aa_gamma_bb_a0');
    ccode(diff(diff(data.functional_a0,gamma_ab,1),gamma_ab,1),'file','v_gamma_ab_gamma_ab_a0');
    ccode(diff(diff(data.functional_a0,gamma_ab,1),gamma_bb,1),'file','v_gamma_ab_gamma_bb_a0');
    ccode(diff(diff(data.functional_a0,gamma_bb,1),gamma_bb,1),'file','v_gamma_bb_gamma_bb_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),tau_a,1),'file','v_rho_a_tau_a_a0');
    ccode(diff(diff(data.functional_a0,rho_a,1),tau_b,1),'file','v_rho_a_tau_b_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),tau_a,1),'file','v_rho_b_tau_a_a0');
    ccode(diff(diff(data.functional_a0,rho_b,1),tau_b,1),'file','v_rho_b_tau_b_a0');
    ccode(diff(diff(data.functional_a0,tau_a,1),tau_a,1),'file','v_tau_a_tau_a_a0');
    ccode(diff(diff(data.functional_a0,tau_a,1),tau_b,1),'file','v_tau_a_tau_b_a0');
    ccode(diff(diff(data.functional_a0,tau_b,1),tau_b,1),'file','v_tau_b_tau_b_a0');
    ccode(diff(diff(data.functional_a0,gamma_aa,1),tau_a,1),'file','v_gamma_aa_tau_a_a0');
    ccode(diff(diff(data.functional_a0,gamma_ab,1),tau_a,1),'file','v_gamma_ab_tau_a_a0');
    ccode(diff(diff(data.functional_a0,gamma_bb,1),tau_a,1),'file','v_gamma_bb_tau_a_a0');
    ccode(diff(diff(data.functional_a0,gamma_aa,1),tau_b,1),'file','v_gamma_aa_tau_b_a0');
    ccode(diff(diff(data.functional_a0,gamma_ab,1),tau_b,1),'file','v_gamma_ab_tau_b_a0');
    ccode(diff(diff(data.functional_a0,gamma_bb,1),tau_b,1),'file','v_gamma_bb_tau_b_a0');

end

%========== ORDER 0 ==========%

try
    ccode(data.functional_b0,'file','functional_b0');
catch exception
    syms qqqq
    ccode(diff(data.functional_b0,qqqq),'file','functional_b0')
end

%========== ORDER 1 ==========%

ccode(diff(data.functional_b0,rho_a,1),'file','v_rho_a_b0');
ccode(diff(data.functional_b0,rho_b,1),'file','v_rho_b_b0');
ccode(diff(data.functional_b0,gamma_aa,1),'file','v_gamma_aa_b0');
ccode(diff(data.functional_b0,gamma_ab,1),'file','v_gamma_ab_b0');
ccode(diff(data.functional_b0,gamma_bb,1),'file','v_gamma_bb_b0');
ccode(diff(data.functional_b0,tau_a,1),'file','v_tau_a_b0');
ccode(diff(data.functional_b0,tau_b,1),'file','v_tau_b_b0');

%========== ORDER 2 ==========%

if (data.deriv2)

    ccode(diff(diff(data.functional_b0,rho_a,1),rho_a,1),'file','v_rho_a_rho_a_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),rho_b,1),'file','v_rho_a_rho_b_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),rho_b,1),'file','v_rho_b_rho_b_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),gamma_aa,1),'file','v_rho_a_gamma_aa_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),gamma_ab,1),'file','v_rho_a_gamma_ab_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),gamma_bb,1),'file','v_rho_a_gamma_bb_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),gamma_aa,1),'file','v_rho_b_gamma_aa_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),gamma_ab,1),'file','v_rho_b_gamma_ab_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),gamma_bb,1),'file','v_rho_b_gamma_bb_b0');
    ccode(diff(diff(data.functional_b0,gamma_aa,1),gamma_aa,1),'file','v_gamma_aa_gamma_aa_b0');
    ccode(diff(diff(data.functional_b0,gamma_aa,1),gamma_ab,1),'file','v_gamma_aa_gamma_ab_b0');
    ccode(diff(diff(data.functional_b0,gamma_aa,1),gamma_bb,1),'file','v_gamma_aa_gamma_bb_b0');
    ccode(diff(diff(data.functional_b0,gamma_ab,1),gamma_ab,1),'file','v_gamma_ab_gamma_ab_b0');
    ccode(diff(diff(data.functional_b0,gamma_ab,1),gamma_bb,1),'file','v_gamma_ab_gamma_bb_b0');
    ccode(diff(diff(data.functional_b0,gamma_bb,1),gamma_bb,1),'file','v_gamma_bb_gamma_bb_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),tau_a,1),'file','v_rho_a_tau_a_b0');
    ccode(diff(diff(data.functional_b0,rho_a,1),tau_b,1),'file','v_rho_a_tau_b_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),tau_a,1),'file','v_rho_b_tau_a_b0');
    ccode(diff(diff(data.functional_b0,rho_b,1),tau_b,1),'file','v_rho_b_tau_b_b0');
    ccode(diff(diff(data.functional_b0,tau_a,1),tau_a,1),'file','v_tau_a_tau_a_b0');
    ccode(diff(diff(data.functional_b0,tau_a,1),tau_b,1),'file','v_tau_a_tau_b_b0');
    ccode(diff(diff(data.functional_b0,tau_b,1),tau_b,1),'file','v_tau_b_tau_b_b0');
    ccode(diff(diff(data.functional_b0,gamma_aa,1),tau_a,1),'file','v_gamma_aa_tau_a_b0');
    ccode(diff(diff(data.functional_b0,gamma_ab,1),tau_a,1),'file','v_gamma_ab_tau_a_b0');
    ccode(diff(diff(data.functional_b0,gamma_bb,1),tau_a,1),'file','v_gamma_bb_tau_a_b0');
    ccode(diff(diff(data.functional_b0,gamma_aa,1),tau_b,1),'file','v_gamma_aa_tau_b_b0');
    ccode(diff(diff(data.functional_b0,gamma_ab,1),tau_b,1),'file','v_gamma_ab_tau_b_b0');
    ccode(diff(diff(data.functional_b0,gamma_bb,1),tau_b,1),'file','v_gamma_bb_tau_b_b0');

end

%========== ORDER 0 ==========%

try
    ccode(data.functional_a0b0,'file','functional_a0b0');
catch exception
    syms qqqq
    ccode(diff(data.functional_a0b0,qqqq),'file','functional_a0b0')
end

%========== ORDER 1 ==========%

ccode(diff(data.functional_a0b0,rho_a,1),'file','v_rho_a_a0b0');
ccode(diff(data.functional_a0b0,rho_b,1),'file','v_rho_b_a0b0');
ccode(diff(data.functional_a0b0,gamma_aa,1),'file','v_gamma_aa_a0b0');
ccode(diff(data.functional_a0b0,gamma_ab,1),'file','v_gamma_ab_a0b0');
ccode(diff(data.functional_a0b0,gamma_bb,1),'file','v_gamma_bb_a0b0');
ccode(diff(data.functional_a0b0,tau_a,1),'file','v_tau_a_a0b0');
ccode(diff(data.functional_a0b0,tau_b,1),'file','v_tau_b_a0b0');

%========== ORDER 2 ==========%

if (data.deriv2)

    ccode(diff(diff(data.functional_a0b0,rho_a,1),rho_a,1),'file','v_rho_a_rho_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),rho_b,1),'file','v_rho_a_rho_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),rho_b,1),'file','v_rho_b_rho_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),gamma_aa,1),'file','v_rho_a_gamma_aa_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),gamma_ab,1),'file','v_rho_a_gamma_ab_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),gamma_bb,1),'file','v_rho_a_gamma_bb_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),gamma_aa,1),'file','v_rho_b_gamma_aa_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),gamma_ab,1),'file','v_rho_b_gamma_ab_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),gamma_bb,1),'file','v_rho_b_gamma_bb_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_aa,1),gamma_aa,1),'file','v_gamma_aa_gamma_aa_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_aa,1),gamma_ab,1),'file','v_gamma_aa_gamma_ab_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_aa,1),gamma_bb,1),'file','v_gamma_aa_gamma_bb_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_ab,1),gamma_ab,1),'file','v_gamma_ab_gamma_ab_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_ab,1),gamma_bb,1),'file','v_gamma_ab_gamma_bb_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_bb,1),gamma_bb,1),'file','v_gamma_bb_gamma_bb_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),tau_a,1),'file','v_rho_a_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_a,1),tau_b,1),'file','v_rho_a_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),tau_a,1),'file','v_rho_b_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,rho_b,1),tau_b,1),'file','v_rho_b_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,tau_a,1),tau_a,1),'file','v_tau_a_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,tau_a,1),tau_b,1),'file','v_tau_a_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,tau_b,1),tau_b,1),'file','v_tau_b_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_aa,1),tau_a,1),'file','v_gamma_aa_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_ab,1),tau_a,1),'file','v_gamma_ab_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_bb,1),tau_a,1),'file','v_gamma_bb_tau_a_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_aa,1),tau_b,1),'file','v_gamma_aa_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_ab,1),tau_b,1),'file','v_gamma_ab_tau_b_a0b0');
    ccode(diff(diff(data.functional_a0b0,gamma_bb,1),tau_b,1),'file','v_gamma_bb_tau_b_a0b0');

end
%========== ORDER 0 ==========%

replaceInFile('UKS_FUNCTIONAL',buildUKS('functional'),[name '_functional.cc']);

%========== ORDER 1 ==========%

replaceInFile('UKS_V1_RHO_A',buildUKS('v_rho_a'),[name '_functional.cc']);
replaceInFile('UKS_V1_RHO_B',buildUKS('v_rho_b'),[name '_functional.cc']);
replaceInFile('UKS_V1_GAMMA_AA',buildUKS('v_gamma_aa'),[name '_functional.cc']);
replaceInFile('UKS_V1_GAMMA_AB',buildUKS('v_gamma_ab'),[name '_functional.cc']);
replaceInFile('UKS_V1_GAMMA_BB',buildUKS('v_gamma_bb'),[name '_functional.cc']);
replaceInFile('UKS_V1_TAU_A',buildUKS('v_tau_a'),[name '_functional.cc']);
replaceInFile('UKS_V1_TAU_B',buildUKS('v_tau_b'),[name '_functional.cc']);

%========== ORDER 2 ==========%

if (data.deriv2)

    replaceInFile('UKS_V2_RHO_A_RHO_A',buildUKS('v_rho_a_rho_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_RHO_B',buildUKS('v_rho_a_rho_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_RHO_B',buildUKS('v_rho_b_rho_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_AA',buildUKS('v_rho_a_gamma_aa'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_AB',buildUKS('v_rho_a_gamma_ab'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_BB',buildUKS('v_rho_a_gamma_bb'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_AA',buildUKS('v_rho_b_gamma_aa'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_AB',buildUKS('v_rho_b_gamma_ab'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_BB',buildUKS('v_rho_b_gamma_bb'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_AA',buildUKS('v_gamma_aa_gamma_aa'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_AB',buildUKS('v_gamma_aa_gamma_ab'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_BB',buildUKS('v_gamma_aa_gamma_bb'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_GAMMA_AB',buildUKS('v_gamma_ab_gamma_ab'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_GAMMA_BB',buildUKS('v_gamma_ab_gamma_bb'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_GAMMA_BB',buildUKS('v_gamma_bb_gamma_bb'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_TAU_A',buildUKS('v_rho_a_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_TAU_B',buildUKS('v_rho_a_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_TAU_A',buildUKS('v_rho_b_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_TAU_B',buildUKS('v_rho_b_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_A_TAU_A',buildUKS('v_tau_a_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_A_TAU_B',buildUKS('v_tau_a_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_B_TAU_B',buildUKS('v_tau_b_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_TAU_A',buildUKS('v_gamma_aa_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_TAU_A',buildUKS('v_gamma_ab_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_TAU_A',buildUKS('v_gamma_bb_tau_a'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_TAU_B',buildUKS('v_gamma_aa_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_TAU_B',buildUKS('v_gamma_ab_tau_b'),[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_TAU_B',buildUKS('v_gamma_bb_tau_b'),[name '_functional.cc']);

else

    replaceInFile('UKS_V2_RHO_A_RHO_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_RHO_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_RHO_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_AA','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_AB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_GAMMA_BB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_AA','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_AB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_GAMMA_BB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_AA','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_AB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_GAMMA_BB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_GAMMA_AB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_GAMMA_BB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_GAMMA_BB','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_A_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_RHO_B_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_A_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_A_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_TAU_B_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_TAU_A','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AA_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_AB_TAU_B','',[name '_functional.cc']);
    replaceInFile('UKS_V2_GAMMA_BB_TAU_B','',[name '_functional.cc']);

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       RKS FUNCTIONALS
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rks_functional = data.functional;
rks_functional_0 = data.functional_a0b0;

%========== ORDER 0 ==========%

% zk <- zk
try
    rksCode(data.functional,'functional');
catch 
    rksCode(simplify(rks_functional),'functional');
end

%========== ORDER 1 ==========%

% v_rho_a <- v_rho_a
rksCode(diff(rks_functional,rho_a,1),'v_rho_a');
% v_gamma_aa <- 2.0 * v_gamma_aa + v_gamma_ab
rksCode(2*diff(rks_functional,gamma_aa,1)+diff(rks_functional,gamma_ab,1),'v_gamma_aa');
% v_tau_a <- v_tau_a
rksCode(diff(rks_functional,tau_a,1),'v_tau_a');

%========== ORDER 2 ==========%

if (data.deriv2)

% v_rho_a_rho_a <- v_rho_a_rho_a + v_rho_a_rho_b 
rksCode(diff(diff(rks_functional,rho_a,1),rho_a,1)+diff(diff(rks_functional,rho_a,1),rho_b,1),'v_rho_a_rho_a');
% v_rho_a_gamma_aa <-  v_rho_a_gamma_aa + v_rho_a_gamma_ab + v_rho_a_gamma_bb 
rksCode(diff(diff(rks_functional,rho_a,1),gamma_aa,1)+diff(diff(rks_functional,rho_a,1),gamma_ab,1)+...
    diff(diff(rks_functional,rho_a,1),gamma_bb,1),'v_rho_a_gamma_aa');
% v_gamma_aa_gamma_aa <- 2 v_gamma_aa_gamma_aa + 4 v_gamma_aa_gamma_ab + 2 v_gamma_aa_gamma_bb + v_gamma_ab_gamma_ab
rksCode(2*diff(diff(rks_functional,gamma_aa,1),gamma_aa,1) + 4*diff(diff(rks_functional,gamma_aa,1),gamma_ab,1) +...
    2*diff(diff(rks_functional,gamma_aa,1),gamma_bb,1) + diff(diff(rks_functional,gamma_ab,1),gamma_ab,1), ...
    'v_gamma_aa_gamma_aa');
% v_rho_a_tau_a <- v_rho_a_tau_a + v_rho_a_tau_b (I think)
rksCode(diff(diff(rks_functional,rho_a,1),tau_a,1) + diff(diff(rks_functional,rho_a,1),tau_b,1),'v_rho_a_tau_a');
% v_tau_a_tau_a <- v_tau_a_tau_a + v_tau_a_tau_b 
rksCode(diff(diff(rks_functional,tau_a,1),tau_a,1)+diff(diff(rks_functional,tau_a,1),tau_b,1),'v_tau_a_tau_a');
% v_tau_a_gamma_aa <-  v_tau_a_gamma_aa + v_tau_a_gamma_ab + v_tau_a_gamma_bb 
rksCode(diff(diff(rks_functional,tau_a,1),gamma_aa,1)+diff(diff(rks_functional,tau_a,1),gamma_ab,1)+...
    diff(diff(rks_functional,tau_a,1),gamma_bb,1),'v_gamma_aa_tau_a');

end

%========== ORDER 0 ==========%

% zk <- zk
rksCode(rks_functional_0,'functional_0');

%========== ORDER 1 ==========%

% v_rho_a <- v_rho_a
rksCode(diff(rks_functional_0,rho_a,1),'v_rho_a_0');
% v_gamma_aa <- 2.0 * v_gamma_aa + v_gamma_ab
rksCode(2*diff(rks_functional_0,gamma_aa,1)+diff(rks_functional_0,gamma_ab,1),'v_gamma_aa_0');
% v_tau_a <- v_tau_a
rksCode(diff(rks_functional_0,tau_a,1),'v_tau_a_0');

%========== ORDER 2 ==========%

if (data.deriv2)

% v_rho_a_rho_a <- v_rho_a_rho_a + v_rho_a_rho_b 
rksCode(diff(diff(rks_functional_0,rho_a,1),rho_a,1)+diff(diff(rks_functional_0,rho_a,1),rho_b,1),'v_rho_a_rho_a_0');
% v_rho_a_gamma_aa <-  v_rho_a_gamma_aa + v_rho_a_gamma_ab + v_rho_a_gamma_bb 
rksCode(diff(diff(rks_functional_0,rho_a,1),gamma_aa,1)+diff(diff(rks_functional_0,rho_a,1),gamma_ab,1)+...
    diff(diff(rks_functional_0,rho_a,1),gamma_bb,1),'v_rho_a_gamma_aa_0');
% v_gamma_aa_gamma_aa <- 2 v_gamma_aa_gamma_aa + 4 v_gamma_aa_gamma_ab + 2 v_gamma_aa_gamma_bb + v_gamma_ab_gamma_ab
rksCode(2*diff(diff(rks_functional_0,gamma_aa,1),gamma_aa,1) + 4*diff(diff(rks_functional_0,gamma_aa,1),gamma_ab,1) +...
    2*diff(diff(rks_functional_0,gamma_aa,1),gamma_bb,1) + diff(diff(rks_functional_0,gamma_ab,1),gamma_ab,1), ...
    'v_gamma_aa_gamma_aa_0');
% v_rho_a_tau_a <- v_rho_a_tau_a + v_rho_a_tau_b (I think)
rksCode(diff(diff(rks_functional_0,rho_a,1),tau_a,1) + diff(diff(rks_functional_0,rho_a,1),tau_b,1),'v_rho_a_tau_a_0');
% v_tau_a_tau_a <- v_tau_a_tau_a + v_tau_a_tau_b 
rksCode(diff(diff(rks_functional_0,tau_a,1),tau_a,1)+diff(diff(rks_functional_0,tau_a,1),tau_b,1),'v_tau_a_tau_a_0');
% v_tau_a_gamma_aa <-  v_tau_a_gamma_aa + v_tau_a_gamma_ab + v_tau_a_gamma_bb 
rksCode(diff(diff(rks_functional_0,tau_a,1),gamma_aa,1)+diff(diff(rks_functional_0,tau_a,1),gamma_ab,1)+...
    diff(diff(rks_functional_0,tau_a,1),gamma_bb,1),'v_gamma_aa_tau_a_0');

end

%========== ORDER 0 ==========%

replaceInFile('RKS_FUNCTIONAL',buildRKS('functional'),[name '_functional.cc']);

%========== ORDER 1 ==========%

replaceInFile('RKS_V1_RHO_A',buildRKS('v_rho_a'),[name '_functional.cc']);
replaceInFile('RKS_V1_GAMMA_AA',buildRKS('v_gamma_aa'),[name '_functional.cc']);
replaceInFile('RKS_V1_TAU_A',buildRKS('v_tau_a'),[name '_functional.cc']);

%========== ORDER 2 ==========%

if (data.deriv2)
    
    replaceInFile('RKS_V2_RHO_A_RHO_A',buildRKS('v_rho_a_rho_a'),[name '_functional.cc']);
    replaceInFile('RKS_V2_RHO_A_GAMMA_AA',buildRKS('v_rho_a_gamma_aa'),[name '_functional.cc']);
    replaceInFile('RKS_V2_GAMMA_AA_GAMMA_AA',buildRKS('v_gamma_aa_gamma_aa'),[name '_functional.cc']);
    replaceInFile('RKS_V2_RHO_A_TAU_A',buildRKS('v_rho_a_tau_a'),[name '_functional.cc']);
    replaceInFile('RKS_V2_TAU_A_TAU_A',buildRKS('v_tau_a_tau_a'),[name '_functional.cc']);
    replaceInFile('RKS_V2_GAMMA_AA_TAU_A',buildRKS('v_gamma_aa_tau_a'),[name '_functional.cc']);

else 

    replaceInFile('RKS_V2_RHO_A_RHO_A','',[name '_functional.cc']);
    replaceInFile('RKS_V2_RHO_A_GAMMA_AA','',[name '_functional.cc']);
    replaceInFile('RKS_V2_GAMMA_AA_GAMMA_AA','',[name '_functional.cc']);
    replaceInFile('RKS_V2_RHO_A_TAU_A','',[name '_functional.cc']);
    replaceInFile('RKS_V2_TAU_A_TAU_A','',[name '_functional.cc']);
    replaceInFile('RKS_V2_GAMMA_AA_TAU_A','',[name '_functional.cc']);

end
