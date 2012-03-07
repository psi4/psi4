function code = codeTree(tree, root, array, expand, trace, use_double)

if (nargin < 6)
    use_double = true;
end
if (nargin < 5)
    trace = false;
end
if (nargin < 4)
    expand = false;
end
if (nargin < 3)
    array = root;
end


newline = sprintf('\n');
code = ['/** Tree Code for ' array ' value. **/' newline];

if (trace)
    code = [code '/**' newline];
    code = [code printTree(tree)];
    code = [code '**/' newline];
end

if (expand)
    tree = expandTree(tree,root);
end

ind = determineDependencies(tree, root);
code = [code newline];

for j = 1:length(ind) - 1
    k = ind(j);
    ccode(tree(k).expr,'file',[tree(k).name '_0']);
    code = [code '// Intermediate ' tree(k).name ':' newline]; 
    code = [code cleanTree([tree(k).name '_0'],tree(k).name, '=', 'double')];
    code = [code newline];
end
ccode(tree(ind(end)).expr,'file',[tree(ind(end)).name '_0']);
code = [code '// Final ' array ':' newline]; 
code = [code cleanTree([tree(ind(end)).name '_0'], array , '=', '')];
code = [code newline];

code = cleanCode(code); 

