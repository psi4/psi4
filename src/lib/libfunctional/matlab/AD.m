function [trees] = AD(tree, root, params, primaries)

if (nargin < 4)
    primaries = {'rho_a' 'rho_b' 'tau_a' 'tau_b' 'gamma_aa' 'gamma_ab' 'gamma_bb'};
end
if (nargin < 3)
    params = {};
end
if (nargin < 2)
    root = 'functional';
end

% Put tree in order
ind = determineDependencies(tree, root);
tree = tree(ind);

% Oth order
ind = strcmp({tree(:).name}, 'functional');
inter0 = tree(find(~ind));
fin0 = tree(find(ind));

old_length = length(tree);
tree = dressTree(tree);
prims = tree(1:(length(tree)-old_length));
primaries_unsorted = { prims(:).name};
primaries_sorted = {}; 
for k = 1:length(primaries)
    if (any(strcmp(primaries{k}, primaries_unsorted)))
        primaries_sorted{end+1} = primaries{k};
    end
end
primaries = primaries_sorted;

% 1st order
old_length = length(tree);
tree = reverseAccumulate(tree,root,params,primaries);

new_part = tree(old_length + 1 : length(tree));
new_names = {new_part(:).name};
ind = zeros(size(new_part));
for k = 1:length(primaries)
    partial = [root '_' primaries{k}];
    ind(find(strcmp(partial,new_names))) = 1;
end

inter1 = new_part(find(~ind));
fin1 = new_part(find(ind)); 

% 2nd order
inter2 = struct([]);
inter2(1).name = 'null';
inter2(1).expr = sym(0);
inter2(1).depend = {};
inter2_names = {};
fin2 = struct([]);
fin2(1).name = 'null';
fin2(1).expr = sym(0);
fin2(1).depend = {};
fin2_names = {};
for v = 1:length(primaries)
    var1 = primaries{v};
    unrolled = removeOrphans(tree, [root '_' var1]);
    old_length = length(unrolled);
    unrolled = reverseAccumulate(unrolled,[root '_' var1], params, primaries);
    new_part = unrolled(old_length + 1 : length(unrolled)); 
    new_names = {new_part(:).name};

    ind = zeros(size(new_part));
    for k = 1:length(primaries)
        partial2 = [root '_' primaries{v} '_' primaries{k}];
        hit = find(strcmp(partial2,new_names));
        ind(hit) = 1;
        if (k < v)
            if (~isempty(hit))
                new_part(hit).name = [root '_' primaries{k} '_' primaries{v}];
                new_names{hit} = [root '_' primaries{k} '_' primaries{v}];
            end
        end
    end
     
    for k = 1:length(ind)
        element = new_part(k);
        if (ind(k) == 1)
            %fin
            if (any(strcmp(fin2_names, element.name))) 
            else
                fin2(end+1) = element;
                fin2_names{end+1} = element.name;
            end
        elseif (ind(k) == 0)
            %inter
            if (any(strcmp(inter2_names, element.name))) 
            else
                inter2(end+1) = element;
                inter2_names{end+1} = element.name;
            end
        end
    end
end

inter2 = inter2(2:end);
fin2 = fin2(2:end);

%disp(sprintf('Primitives:'))
%printTree(prims);
%disp(sprintf('0th order intermediates:'))
%printTree(inter0);
%disp(sprintf('1st order intermediates:'))
%printTree(inter1);
%disp(sprintf('2nd order intermediates:'))
%printTree(inter2);
%disp(sprintf('0th order final quantities:'))
%printTree(fin0);
%disp(sprintf('1st order final quantities:'))
%printTree(fin1);
%disp(sprintf('2nd order final quantities:'))
%printTree(fin2);

trees.prims = prims;
trees.inter0 = inter0;
trees.inter1 = inter1;
trees.inter2 = inter2;
trees.fin0 = fin0;
trees.fin1 = fin1;
trees.fin2 = fin2;

