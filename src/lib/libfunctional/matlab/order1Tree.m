function tree2 = order1Tree(tree, root, var1, params, primaries)

valid1 = false;
for k = 1:length(tree)
    var_str = findsym(tree(k).expr);    
    vars = {};
    [T var_str] = strtok(var_str, ',');
    while (~isempty(T))
        vars{end+1} = T;
        [T var_str] = strtok(var_str, ',');
    end
    if (any(strcmp(vars, var1)))    
        valid1 = true;
        break;
    end
end
if (~valid1)
    tree2(1).name = [root '_' var1];
    tree2(1).depend = {};
    tree2(1).expr = sym('0');
    return;
end

if (nargin < 5)
    primaries = {'rho_a' 'rho_b' 'tau_a' 'tau_b' 'gamma_aa' 'gamma_bb' 'gamma_ab'};
end
if (nargin < 4)
    params = {};
end

tree = removeOrphans(tree,root);
ind = determineDependencies(tree, root);

tree2 = struct([]);

for k = 1:length(ind) - 1
    
    index = ind(k);
    name = tree(index).name;

    % value 0
    tree2(end+1).name = name;
    tree2(end).expr = tree(index).expr + 0;

    for l = 1:length(tree(index).depend)
        depend = tree(index).depend{l};

        tree2(end+1).name = [name '_' depend];
        tree2(end).expr = diff(tree(index).expr,sym(depend));

    end
        
    tree2(end+1).name = [name '_' var1];
    tree2(end).expr = diff(tree(index).expr,sym(var1));
        
    for l = 1:length(tree(index).depend)
        depend = tree(index).depend{l};
        
        tree2(end).expr = tree2(end).expr + sym([name '_' depend]) * sym([depend '_' var1]);

    end

end

index = ind(end);
name = tree(index).name;

for l = 1:length(tree(index).depend)
    depend = tree(index).depend{l};

    tree2(end+1).name = [name '_' depend];
    tree2(end).expr = diff(tree(index).expr,sym(depend));

end
    
tree2(end+1).name = [name '_' var1];
tree2(end).expr = diff(tree(index).expr,sym(var1));
    
for l = 1:length(tree(index).depend)
    depend = tree(index).depend{l};
    
    tree2(end).expr = tree2(end).expr + sym([name '_' depend]) * sym([depend '_' var1]);

end

tree2 = buildDependencies(tree2,params,primaries);
tree2 = removeNumbers(tree2);
tree2 = buildDependencies(tree2,params,primaries);
tree2 = removeOrphans(tree2,[name '_' var1]);

