function s = buildRKSII(code)

s = '';

newline = sprintf('\n');
code = [newline code];
code = regexprep(code,'\n','\n    ');
code = code(2:end-4);

s = [s 'if (rho_a[index] > cutoff_) {' newline];
s = [s code];
s = [s '} ' newline];

