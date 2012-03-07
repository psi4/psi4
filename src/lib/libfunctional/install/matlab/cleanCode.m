function code = cleanCode(code)

code = regexprep(code,'([^_])rho_a','$1rho_a[index]');
code = regexprep(code,'([^_])rho_b','$1rho_b[index]');
code = regexprep(code,'([^_])gamma_aa','$1gamma_aa[index]');
code = regexprep(code,'([^_])gamma_ab','$1gamma_ab[index]');
code = regexprep(code,'([^_])gamma_bb','$1gamma_bb[index]');
code = regexprep(code,'([^_])tau_a','$1tau_a[index]');
code = regexprep(code,'([^_])tau_b','$1tau_b[index]');

