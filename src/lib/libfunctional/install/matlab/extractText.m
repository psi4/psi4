function ret = extractText(filename)

%Get text (including newlines) and 
% discard file
file = fopen(filename,'r');

ret = '';
token = fgets(file);
while ischar(token)
    ret = [ret token];
    token = fgets(file);
end

fclose(file);

rm_str = sprintf('rm %s',filename);
system(rm_str);
