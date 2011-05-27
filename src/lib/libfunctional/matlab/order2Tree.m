function tree2 = order2Tree(tree, root, var1, var2, params, primaries)

valid1 = false;
valid2 = false;
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
for k = 1:length(tree)
    var_str = findsym(tree(k).expr);    
    vars = {};
    [T var_str] = strtok(var_str, ',');
    while (~isempty(T))
        vars{end+1} = T;
        [T var_str] = strtok(var_str, ',');
    end
    if (any(strcmp(vars, var2)))    
        valid2 = true;
        break;
    end
end
if (~valid1 || ~ valid2)
    tree2(1).name = [root '_' var1 '_' var2];
    tree2(1).depend = {};
    tree2(1).expr = sym('0');
    return;
end

if (nargin < 6)
    primaries = {'rho_a' 'rho_b' 'tau_a' 'tau_b' 'gamma_aa' 'gamma_bb' 'gamma_ab'};
end
if (nargin < 5)
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

    % value 1
    for l = 1:length(tree(index).depend)
        depend1 = tree(index).depend{l};

        tree2(end+1).name = [name '_' depend1];
        tree2(end).expr = diff(tree(index).expr,sym(depend1));

    end
        
    tree2(end+1).name = [name '_' var1];
    tree2(end).expr = diff(tree(index).expr,sym(var1));
    for l = 1:length(tree(index).depend)
        depend1 = tree(index).depend{l};
        
        tree2(end).expr = tree2(end).expr + sym([name '_' depend1]) * sym([depend1 '_' var1]);

    end

    if (~strcmp(var1,var2))
        tree2(end+1).name = [name '_' var2];
        tree2(end).expr = diff(tree(index).expr,sym(var2));
        for l = 1:length(tree(index).depend)
            depend1 = tree(index).depend{l};
            
            tree2(end).expr = tree2(end).expr + sym([name '_' depend1]) * sym([depend1 '_' var2]);
    
        end
    end

    % value 2
    for l = 1:length(tree(index).depend)
        depend1 = tree(index).depend{l};

        for m = 1:l
            depend2 = tree(index).depend{m};

            tree2(end+1).name = [name '_' depend1 '_' depend2];
            tree2(end).expr = diff(diff(tree(index).expr,sym(depend1)),sym(depend2));
        end
    end

    tree2(end+1).name = [name '_' var1 '_' var2];
    tree2(end).expr = diff(diff(tree(index).expr,sym(var1)),sym(var2));

    for l = 1:length(tree(index).depend)
        depend1 = tree(index).depend{l};

        tree2(end).expr = tree2(end).expr + sym([name '_' depend1]) * sym([depend1 '_' var1 '_' var2]);
    end
    
    for l = 1:length(tree(index).depend)
        depend1 = tree(index).depend{l};

        for m = 1:l
            depend2 = tree(index).depend{m};

            tree2(end).expr = tree2(end).expr + sym([name '_' depend1 '_' depend2]) * sym([depend1 '_' var1]) * sym([depend2 '_' var2]);
            if (~strcmp(depend1,depend2))
                tree2(end).expr = tree2(end).expr + sym([name '_' depend1 '_' depend2]) * sym([depend2 '_' var1]) * sym([depend1 '_' var2]);
            end
        end
    end

end

index = ind(end);
name = tree(index).name;

% value 2
for l = 1:length(tree(index).depend)
    depend1 = tree(index).depend{l};

    tree2(end+1).name = [name '_' depend1];
    tree2(end).expr = diff(tree(index).expr,sym(depend1));

end

for l = 1:length(tree(index).depend)
    depend1 = tree(index).depend{l};

    for m = 1:l
        depend2 = tree(index).depend{m};

        tree2(end+1).name = [name '_' depend1 '_' depend2];
        tree2(end).expr = diff(diff(tree(index).expr,sym(depend1)),sym(depend2));
    end
end

tree2(end+1).name = [name '_' var1 '_' var2];
tree2(end).expr = diff(diff(tree(index).expr,sym(var1)),sym(var2));

for l = 1:length(tree(index).depend)
    depend1 = tree(index).depend{l};

    tree2(end).expr = tree2(end).expr + sym([name '_' depend1]) * sym([depend1 '_' var1 '_' var2]);
end

for l = 1:length(tree(index).depend)
    depend1 = tree(index).depend{l};

    for m = 1:l
        depend2 = tree(index).depend{m};

        tree2(end).expr = tree2(end).expr + sym([name '_' depend1 '_' depend2]) * sym([depend1 '_' var1]) * sym([depend2 '_' var2]);
        if (~strcmp(depend1,depend2))
            tree2(end).expr = tree2(end).expr + sym([name '_' depend1 '_' depend2]) * sym([depend2 '_' var1]) * sym([depend1 '_' var2]);
        end
    end
end

tree2 = buildDependencies(tree2,params,primaries);
tree2 = removeNumbers(tree2);
tree2 = buildDependencies(tree2,params,primaries);
tree2 = removeOrphans(tree2,[name '_' var1 '_' var2]);

