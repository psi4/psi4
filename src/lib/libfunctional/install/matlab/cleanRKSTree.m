function tree = cleanRKSTree(tree, root, params)

tree = subsTreeVariable(tree,sym('rho_b'),sym('rho_a'));
tree = subsTreeVariable(tree,sym('gamma_ab'),sym('gamma_aa'));
tree = subsTreeVariable(tree,sym('gamma_bb'),sym('gamma_aa'));
tree = subsTreeVariable(tree,sym('tau_b'),sym('tau_a'));

tree = removeDuplicates(tree);
tree = buildDependencies(tree, params);
tree = removeOrphans(tree,root);
