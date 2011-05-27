function tree1 = reverseAccumulate(tree, root, params, primaries) 

if (nargin < 4)
    primaries = {'rho_a' 'rho_b' 'tau_a' 'tau_b' 'gamma_aa' 'gamma_bb' 'gamma_ab'};
end
if (nargin < 3)
    params = {};
end

tree = removeOrphans(tree,root);
ind = determineDependencies(tree, root);
tree = tree(ind);
ind = determineDependencies(tree, root);
rev_ind = ind(end:-1:1);

parents = struct([]);
for Q = 1:length(tree)
    parents(Q).vals = [];
    name = tree(Q).name;
    for P = 1:length(tree)
        if (any(strcmpi(name, tree(P).depend)))
            parents(Q).vals(end+1) = P;
        end
    end
end

tree1 = struct([]);

% Implicit Seed
%tree1(1).name = [root '_' root];
%tree1(1).expr = sym(1);
%tree1(1).depend = {};

for Qtemp = 1:length(rev_ind)
    Q = rev_ind(Qtemp);
    tree1(Qtemp).name = [root '_' tree(Q).name];
    tree1(Qtemp).expr = 0;
    tree1(Qtemp).depend = {};
    name = tree(Q).name;
    for P = 1:length(parents(Q).vals)
        parent_index = parents(Q).vals(P);
        if (length(parents(parent_index).vals) == 0) 
            tree1(Qtemp).expr = tree1(Qtemp).expr + diff(tree(parent_index).expr, sym(name));
        else
            tree1(Qtemp).expr = tree1(Qtemp).expr + diff(tree(parent_index).expr, sym(name)) * sym([root '_' tree(parent_index).name]);
        end
        if (strcmp(class(tree(Qtemp).expr) ,'double'))
            tree1(Qtemp).expr = sym(num2str(tree(Qtemp).expr));
        end
    end
end

tree1 = tree1(2:end);

for Q = 1:length(tree1)
    var_str = findsym(tree1(Q).expr);    
    vars = {};
    [T var_str] = strtok(var_str, ',');
    while (~isempty(T))
        vars{end+1} = T;
        [T var_str] = strtok(var_str, ',');
    end

    for d = 1:length(vars)
        is_param = false;
        if (any(strcmpi(vars{d}, params)))
            is_param = true;
        end 
        if (strcmpi(vars{d}, tree1(Q).name))
            is_param = true;
        end
        if (~is_param)
            tree1(Q).depend{end+1} = vars{d};
        end
    end
end

tree1 = [tree tree1];
