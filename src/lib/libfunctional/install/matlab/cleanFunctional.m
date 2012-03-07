function polished = cleanFunctional(name, handle)

%Now deal with leading variables
%Newline is a hack
raw = [sprintf('\n') extractText(name)];
raw = regexprep(raw,'rho_a','rho_a[index]');
raw = regexprep(raw,'rho_b','rho_b[index]');
raw = regexprep(raw,'gamma_aa','gamma_aa[index]');
raw = regexprep(raw,'gamma_ab','gamma_ab[index]');
raw = regexprep(raw,'gamma_bb','gamma_bb[index]');
raw = regexprep(raw,'tau_a','tau_a[index]');
raw = regexprep(raw,'tau_b','tau_b[index]');
raw = extractHeavisides(raw);
raw = extractDiracs(raw);

[starts ends extents match tokens] = regexp(raw,'\n\s*?(t)(\d+)');

cleaned = '';
stop = 1;
for k = 1:length(starts)-1
    cleaned = [cleaned raw(stop:extents{k}(1,1)-1) 'double t' tokens{k}{2}];
    stop = extents{k}(2,2) + 1;
end 

k = length(starts);
cleaned = [cleaned raw(stop:extents{k}(1,1)-1) handle  '_[index]' raw(extents{k}(2,2)+1:end)];

newline = sprintf('\n');
cleaned = regexprep(cleaned,'\n  ','\n    ');

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
