function polished = cleanTree(name, handle, operator, prefix)

if (nargin < 3)
    operator = '=';
end
if (nargin < 4)
    prefix = 'double';
end

if (~isempty(prefix))
    prefix = [prefix ' '];
end

%Now deal with leading variables
%Newline is a hack
raw = [sprintf('\n') extractText(name)];

[starts ends extents match tokens] = regexp(raw,'\n\s*?(t)(\d+)');

cleaned = '';
stop = 1;
for k = 1:length(starts)-1
    cleaned = [cleaned raw(stop:extents{k}(1,1)-1) 'double t' tokens{k}{2}];
    stop = extents{k}(2,2) + 1;
end 

k = length(starts);
cleaned = [cleaned raw(stop:extents{k}(1,1)-1) prefix handle ' ' operator raw(extents{k}(2,2)+3:end)];

newline = sprintf('\n');
cleaned = regexprep(cleaned,'\n  ','\n');

%undo Hack
cleaned = cleaned(2:end);

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
    if (delta > 100)
        if (~isempty(strfind(')+*-/ ',cleaned(k))))
            polished = [polished cleaned(left:k) ' \' newline '       '];
            delta = 3;
            left = k + 1;
        end
    end
end

polished = [polished cleaned(left:end)];
