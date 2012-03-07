function tree = buildDependencies(tree, params, primaries)

if (nargin < 3)
    primaries = {'rho_a' 'rho_b' 'tau_a' 'tau_b' 'gamma_aa' 'gamma_bb' 'gamma_ab'};
end
if (nargin < 2)
    params = {};
end

reserved = { params{:} primaries{:} };
names = {tree(:).name};

if (isfield(tree(1),'depend'))
    tree = rmfield(tree,'depend');
end
for k = 1:length(tree)
    tree(k).depend = {};
end

for k = 1:length(tree)
    string = findsym(tree(k).expr);

    symvars = {};
    counter = 1;    

    [T string] = strtok(string,',');
    while (~isempty(T)) 
        symvars{counter} = T;
        counter = counter + 1;
        [T string] = strtok(string,',');
    end

    for l = 1:length(symvars)
        var = symvars{l};
        if (~any(strcmp(var,reserved)))
            % Its a depend
            if (~any(strcmp(var,names)))
                %And its not in the tree
                error(['Dependent variable ' var ' is required by ' tree(k).name ', but is not in the tree.']);
            end 
            tree(k).depend{end+1} = var;   
        end
    end
    
end
