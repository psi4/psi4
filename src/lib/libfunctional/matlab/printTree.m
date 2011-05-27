function s2 = printTree(tree)

s = sprintf('Function Tree: %d elements.\n',length(tree)) ;
newline = sprintf('\n');
s = [ s '------------------------' newline];

for k = 1:length(tree)
    dep = '';
    for l = 1:length(tree(k).depend)
        dep = [dep tree(k).depend{l} ','];
    end

    if (isempty(dep))
        dep = 'None';
    else
        dep = dep(1:end-1);
    end

    s = [ s 'Element:      ' sprintf('%d',k) newline];
    s = [ s 'Name:         ' tree(k).name newline];
    cleaned = char(tree(k).expr);    
    %Clean up lines (wrap at 100 cols, indenting by 3)
    polished = '';
    
    left = 1;
    delta = 0;
    for k = 1:length(cleaned)
        if (cleaned(k) == newline)
            delta = 0;
            continue;
        else
            delta = delta + 1;
        end
        if (delta > 80)
            if (~isempty(strfind(')+*-/ ',cleaned(k))))
                polished = [polished cleaned(left:k) ' ...' newline '       '];
                delta = 3;
                left = k + 1;
            end
        end
    end
    
    expr = [polished cleaned(left:end)];

    s = [ s 'Expression:   ' expr newline];
    s = [ s 'Dependencies: ' dep newline];
    s = [ s '------------------------' newline];
    
end

if (nargout < 1)
    disp(s)
else
    s2 = s;
end
