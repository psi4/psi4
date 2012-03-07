function indices = resolveDependencies(indices)

% Remove redundancies
% Keep 'em all at first
good_list = ones(size(indices));

% Check for early repeats
for k = length(indices):-1:1
    val = indices(k);
    good_list(find(indices(1:k-1) == val)) = 0;
end

indices = indices(find(good_list));
