function indices = determineDependencies(tree, root)

root_index = find(strcmp({tree(:).name}, root));

% Find all dependencies
indices = visitDependencies(tree,[],root_index);

% Remove redundant dependencies
indices = resolveDependencies(indices);

% Reverse order to put in priority order
indices = indices(end:-1:1);
