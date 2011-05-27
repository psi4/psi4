function tree2 = addTrees(trees, roots, root, weights)

if (nargin < 4)
    weights = ones(length(trees),1);
end

tree2 = struct([]);

depend = {};
for k = 1:length(trees)
    index = find(strcmp({trees{k}(:).name},roots{k}));
    
    for l = 1:length(trees{k}(index).depend)
        depend{end+1} = trees{k}(index).depend{l};
    end

    tree2 = [tree2 trees{k}(1:index-1) trees{k}(index+1:end) ];

end
depend = unique(depend);

tree2(end+1).name = root;
tree2(end).expr = 0;
tree2(end).depend = depend;

for k = 1:length(trees)
    index = find(strcmp({trees{k}(:).name},roots{k}));
    tree2(end).expr = tree2(end).expr + weights(k) * trees{k}(index).expr;
    
end

tree2 = removeDuplicates(tree2);
