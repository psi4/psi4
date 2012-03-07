function tree = removeOrphans(tree, root)
% Could make the tree more visually appealing
% Unfortunately kills the root
while (true)

kill = [];
for k = 1:length(tree)
    key = tree(k).name;
    if (strcmp(root,key))
        continue;
    end
    orphan = true;
    for l = 1:length(tree)
        if (any(strcmp(tree(l).depend,key)))
            orphan = false;
            break;
        end
    end
    if (orphan)
        kill = [kill k];
    end
end

if (isempty(kill))
    break;
end

tree(kill) = [];

end

end

