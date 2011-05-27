function indices = visitDependencies(tree,indices,root_index)

indices = [indices root_index];

for k = 1:length(tree(root_index).depend)
    depend_index = find(strcmp({tree(:).name}, tree(root_index).depend{k}));
    indices = visitDependencies(tree,indices,depend_index);
end
