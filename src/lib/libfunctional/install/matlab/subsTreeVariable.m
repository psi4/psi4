function tree = subsTreeVariable(tree,old,new)

key = findsym(old);
for k = 1:length(tree)
    tree(k).expr = subs(tree(k).expr,old,new,0);
    if (strcmp(class(tree(k).expr),'double'))
        tree(k).expr = sym(num2str(tree(k).expr));
    end
end

tree = removeNumbers(tree);
