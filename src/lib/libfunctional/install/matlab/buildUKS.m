function s = buildUKS(root)

root_a0 = [root '_a0'];
root_b0 = [root '_b0'];
root_a0b0 = [root '_a0b0'];

s = '';
newline = sprintf('\n');
s = [s 'if (rho_a[index] > cutoff_ && rho_b[index] > cutoff_) {' newline];
s = [s cleanFunctional(root, root)];
s = [s '} else if (rho_a[index] > cutoff_) {' newline];
s = [s cleanFunctional(root_b0, root)];
s = [s '} else if (rho_b[index] > cutoff_) {' newline];
s = [s cleanFunctional(root_a0, root)];
s = [s '} else {' newline];
s = [s cleanFunctional(root_a0b0, root)];
s = [s '} ' newline];
