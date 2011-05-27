function s = buildUKSII(code)

newline = sprintf('\n');
code{1} = [newline code{1}];
code{1} = regexprep(code{1},'\n','\n    ');
code{1} = code{1}(2:end-4);
code{2} = [newline code{2}];
code{2} = regexprep(code{2},'\n','\n    ');
code{2} = code{2}(2:end-4);
code{3} = [newline code{3}];
code{3} = regexprep(code{3},'\n','\n    ');
code{3} = code{3}(2:end-4);

s = '';
s = [s 'if (rho_a[index] > cutoff_ && rho_b[index] > cutoff_) {' newline];
s = [s code{1}];
s = [s '} else if (rho_a[index] > cutoff_) {' newline];
s = [s code{3}];
s = [s '} else if (rho_b[index] > cutoff_) {' newline];
s = [s code{2}];
s = [s '} ' newline];
