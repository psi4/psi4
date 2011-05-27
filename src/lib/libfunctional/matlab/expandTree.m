function tree = expandTree(tree, root)

ind = determineDependencies(tree, root);

for k = 1:length(ind) - 1
    index = ind(k);
    name = tree(index).name;
    depends = [];
    for l = 1:length(tree)
        if (any(strcmp(tree(l).depend,name)))
            depends = [depends l];
        end
    end

    for m = 1:length(depends)
        depend = depends(m);
        tree(depend).expr = subs(tree(depend).expr,sym(name), tree(index).expr);
        if (strcmp(class(tree(depend).expr),'double'))
            tree(depend).expr = sym(num2str(tree(depend).expr));
        end
    end
end

tree = tree(ind(end));
tree(1).depend = {};
