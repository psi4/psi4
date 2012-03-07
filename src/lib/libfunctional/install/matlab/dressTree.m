function tree = dressTree(tree, variables)

if (nargin < 2)
    variables = {};
    deps = {};
    for Q = 1:length(tree)
        var_str = findsym(tree(Q).expr);    
        vars = {};
        [T var_str] = strtok(var_str, ',');
        while (~isempty(T))
            vars{end+1} = T;
            [T var_str] = strtok(var_str, ',');
        end
        if(any(strcmpi(vars,'rho_a')))
            variables{end + 1} = 'rho_a';     
        end
        if(any(strcmpi(vars,'rho_b')))
            variables{end + 1} = 'rho_b';     
        end
        if(any(strcmpi(vars,'gamma_aa')))
            variables{end + 1} = 'gamma_aa';     
        end
        if(any(strcmpi(vars,'gamma_ab')))
            variables{end + 1} = 'gamma_ab';     
        end
        if(any(strcmpi(vars,'gamma_bb')))
            variables{end + 1} = 'gamma_bb';     
        end
        if(any(strcmpi(vars,'tau_a')))
            variables{end + 1} = 'tau_a';     
        end
        if(any(strcmpi(vars,'tau_b')))
            variables{end + 1} = 'tau_b';     
        end
    
    end

    variables = unique(variables);

end

prims = struct([]);
for Q = 1:length(variables)
    prims(Q).name = variables{Q};
    prims(Q).expr = sym(variables{Q}); 
    prims(Q).depend = {};
end

tree = [prims tree];

for Q = length(variables) + 1: length(tree)
    var_str = findsym(tree(Q).expr);    
    vars = {};
    [T var_str] = strtok(var_str, ',');
    while (~isempty(T))
        vars{end+1} = T;
        [T var_str] = strtok(var_str, ',');
    end

    for d = 1:length(variables)
        if (any(strcmpi(vars, variables{d})))
            tree(Q).depend{end+1} = variables{d};
        end 
    end
        
end

