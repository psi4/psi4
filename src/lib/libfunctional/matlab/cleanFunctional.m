function polished = cleanFunctional(name, handle)

%Clean internal variables
rho_a_str = sprintf('sed -i s/%s/%s/g "%s"','rho_a','rho_a\[index\]',name);
rho_b_str = sprintf('sed -i s/%s/%s/g "%s"','rho_b','rho_b\[index\]',name);
gamma_aa_str = sprintf('sed -i s/%s/%s/g "%s"','gamma_aa','gamma_aa\[index\]',name);
gamma_ab_str = sprintf('sed -i s/%s/%s/g "%s"','gamma_ab','gamma_ab\[index\]',name);
gamma_bb_str = sprintf('sed -i s/%s/%s/g "%s"','gamma_bb','gamma_bb\[index\]',name);
tau_a_str = sprintf('sed -i s/%s/%s/g "%s"','tau_a','tau_a\[index\]',name);
tau_b_str = sprintf('sed -i s/%s/%s/g "%s"','tau_b','tau_b\[index\]',name);

system(rho_a_str);
system(rho_b_str);
system(gamma_aa_str);
system(gamma_ab_str);
system(gamma_bb_str);
system(tau_a_str);
system(tau_b_str);

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
