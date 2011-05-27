function s = buildRKS(root)

root_0 = [root '_0'];

s = '';
newline = sprintf('\n');
s = [s 'if (rho_a[index] > cutoff_) {' newline];
s = [s cleanFunctional(root, root)];
s = [s '} else {' newline];
s = [s cleanFunctional(root_0, root)];
s = [s '} ' newline];

