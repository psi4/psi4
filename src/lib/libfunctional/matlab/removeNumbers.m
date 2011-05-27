function tree = removeNumbers(tree)

constants = {};
values = [];

for k = 1:length(tree)
    if (isempty(findsym(tree(k).expr)))
        constants{end+1} = tree(k).name;
        values(end+1) = double(tree(k).expr);
    end
end

for l = 1:length(constants)
    tree = subsTreeElement(tree,sym(constants{l}),values(l));
end
