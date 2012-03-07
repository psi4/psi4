function tree = subsTreeElement(tree, old, new, recursive)

if (nargin < 4)
    recursive = false;
end

% Tree is a DAG structure array
% old is a sym 
% new is a sym or a double

key = findsym(old);
index = find(strcmp({tree(:).name}, key));
if (isempty(index))
    % nothing to do
    return;
elseif (length(index) > 1)
    error(['Duplicate element ' key ' in tree.']);
end

% If the replacement is a sym that is an intermediate,
% only that element of the tree is affected. 
% Dependencies are NOT affected and MUST be 
% resolved directly 
if (strcmp(class(new), 'sym') && ~isempty(index))
    tree(index).expr = new;  
    return;
end

% If the replacement is a double, the expression 
% is removed from the tree and embedded in the 
% parent expressions directly. The dependencies are 
% resolved. If this forces any parents to a double, 
% the process continues. If no parents are identified,
% the element is retained as a symbolic constant
if (strcmp(class(new) , 'double'))
    
    % First find the parents
    parents = [];
    for k = 1:length(tree)
        if (any(strcmp(tree(k).depend, key)))
            parents = [parents k];
        end
    end 
    parents = unique(parents); 

    % If there are no parents, retain the entry 
    % as a sym (it could be root), but kill its 
    % dependencies to prevent redundant work in 
    % ccode
    %
    % Otherwise, mark the entry, and substitute
    % the double into the parent expressions
    if (isempty(parents))
        tree(index).expr = sym(num2str(new)); 
        tree(index).depend = {};
    else
        %Substitute the double into the parents after
        %marking the entry for removal
       
        tree(index).expr = 'kill'; 

        % Embed the 
        for k = 1:length(parents)
            parent = parents(k);
            tree(parent).expr = subs(tree(parent).expr, old, new);
            
            % if the subs results in a sym, that's dandy
            % else it too is double, and we repeat the procedure
            if (strcmp(class(tree(parent).expr),'sym') && ~isempty(findsym(tree(parent).expr)))
                tree(parent).depend(find(strcmp(tree(parent).depend, key))) = [];
            else
                tree = subsTreeElement(tree,sym(tree(parent).name), double(tree(parent).expr), true); 
            end
        end
    end

    % The base caller now destroys all marked entries
    if (~recursive)
        kill = [];
        for k = 1:length(tree)
            if (strcmp(class(tree(k).expr), 'char'))
                kill = [kill k];
            elseif (strcmp(class(tree(k).expr), 'sym'))
                if (strcmp(findsym(tree(k).expr), 'kill'))
                    kill = [kill k];
                end
            end
        end
        tree(kill) = [];
    end

    return;
end
end

