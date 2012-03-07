function S = getSuperFunctionals(xlsfile)

[n t r] = xlsread(xlsfile);

headers = {r{1,:}};
r = r(2:end,:);

[rows cols] = size(r);

for k = 1:rows
    for l = 1:cols
        header = headers{l};
        if (strcmp(header,'weights'))
            S(k).(header) = str2num(r{k,l});
        elseif (strcmp(header,'functionals'))
            [t rem] = strtok(r{k,l});
            index = 1;
            funcs{index} = t;
            while (~isempty(rem))
                [t rem] = strtok(rem);
                index = index + 1;
                funcs{index} = t;
            end
            S(k).(header) = funcs;
        else
            S(k).(header) = r{k,l};
        end
    end
end

