function tree = removeDuplicates(tree)

worked = false;

while true

    kill = [];

    for k = 1:length(tree)-1
        for l = k+1:length(tree)
            if (tree(k).expr == tree(l).expr)
                key_a = tree(k).name;
                key_b = tree(l).name;
    
                parents_b = [];
                for m = 1:length(tree)
                    if (any(strcmp(tree(m).depend, key_b)))
                        parents_b = [parents_b m];
                    end
                end
                parents_b = unique(parents_b); 
                
                for m = 1:length(parents_b)
                    parent = parents_b(m);
                    tree(parent).expr = subs(tree(parent).expr,sym(key_b),sym(key_a),0);
                    if (strcmp(class(tree(parent).expr), 'double'))
                        tree(parent).expr = sym(num2str(tree(parent).expr));
                    end
    
                    tree(parent).depend(find(strcmp(tree(parent).depend,key_b))) = [];
                    if (~any(strcmp(tree(parent).depend,key_a)))
                        tree(parent).depend{end+1} = key_a;
                    end
                end
    
                kill = [kill l];
            end 
        end
    end

    if (isempty(kill))
        break;
    end

    worked = true;

    tree(kill) = [];

end

if (worked)
    tree = removeNumbers(tree);
end
