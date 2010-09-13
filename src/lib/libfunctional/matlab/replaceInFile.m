function replaceInFile(old,new,filename)

file = fopen(filename,'r');
temp = fopen('temp','w');

line = fgets(file);

regex = ['(\s*)' old];

while ischar(line)
    [starts ends extents matches tokens] = regexp(line,regex);
    if (isempty(starts))
        fprintf(temp,'%s',line);
    else
        str = tokens{1}{1};
        str = [str regexprep(new,'\n',['\n', tokens{1}{1}])];
        index = length(str);
        for k = length(str):-1:1
            if (str(k) ~= ' ')
                index = k;
                 break;
            end
        end
        str = str(1:index);
        fprintf(temp,'%s',str);
    end
    line = fgets(file);
end

fclose(temp);
fclose(file);

rm_str = sprintf('rm %s',filename);
system(rm_str);

cp_str = sprintf('mv temp %s',filename);
system(cp_str);
