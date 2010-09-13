function createFunctionalLibrary(control,alias,templates,target)

%Move tempaltes and existing files over
mkdir(target)
chdir(templates)

template_str = sprintf('cp *.template ../%s', target);
cc_str = sprintf('cp *.cc ../%s', target);
h_str = sprintf('cp *.h ../%s', target);
in_str = sprintf('cp *.in ../%s', target);
%py_str = sprintf('cp *.py ../%s', target);

system(template_str);
system(cc_str);
system(h_str);
system(in_str);
%system(py_str);

chdir('../')
chdir(target)

% Copy factory template to desired factory.cc file
current_dir = dir('./');

if (any(strcmpi(['factory.cc'],{current_dir.name})))
    rm_str = sprintf('rm "factory.cc"');
    system(rm_str);
end

copy_str = sprintf('cp factory.cc.template factory.cc');
system(copy_str);

%=========== BEGINNING OF FUNCTIONAL ADDITION ===============%

header_str = '';
constructor_str = '';
printer_str_x = '';
printer_str_c = '';

disp(sprintf('Creating Functional Library in directory %s',target));

ctrl = fopen(control);

name = fgetl(ctrl);
while ischar(name)
    disp(sprintf('  Creating Functional %s',name));
    getF = inline(sprintf('get%sFunctional()',name));
    data = getF(0);
    
    createFunctional(data);
    header_str = [header_str sprintf('#include "%s_functional.h"\n',name)];
    constructor_str = [constructor_str sprintf('if (boost::to_upper(name) == "%s")\n',name)];
    constructor_str = [constructor_str sprintf('    return boost::shared_ptr<Functional> (new %s_Functional(npoints,deriv));\n',name)];
    if (data.is_exchange)
        if (data.is_lsda)
            lsda = 'X';
        else
            lsda = ' ';
        end
        if (data.is_gga)
            gga = 'X';
        else
            gga = ' ';
        end
        if (data.is_meta)
            meta = 'X';
        else
            meta = ' ';
        end
        printer_str_x = [printer_str_x sprintf('f << "   %-10s      %s        %s        %s " << endl;\n',name,lsda,gga,meta)];
    else
        if (data.is_lsda)
            lsda = 'X';
        else
            lsda = ' ';
        end
        if (data.is_gga)
            gga = 'X';
        else
            gga = ' ';
        end
        if (data.is_meta)
            meta = 'X';
        else
            meta = ' ';
        end
        printer_str_c = [printer_str_c sprintf('f << "   %-10s      %s        %s        %s " << endl;\n',name,lsda,gga,meta)];
    end
    name = fgetl(ctrl);
end

fclose(ctrl);

replaceInFile('HEADERS',header_str,'factory.cc');
replaceInFile('CONSTRUCTORS',constructor_str,'factory.cc');
replaceInFile('PRINTERS_X',printer_str_x,'factory.cc');
replaceInFile('PRINTERS_C',printer_str_c,'factory.cc');

%============= END OF FUNCTIONAL ADDITION ===============%

%=========== BEGINNING OF SUPERFUNCTIONAL ADDITION ===============%

alias_ctrl = fopen(alias);
alias_str = '';

name = fgetl(alias_ctrl);
while ischar(name)
    disp(sprintf('  Creating SuperFunctional %s',name));
    getF = inline(sprintf('get%sSuperFunctional()',name));
    data = getF(0);

    dashD_treatment = 'NULL';
    if (~isempty(data.dashD))
        dashD_treatment = sprintf('Dispersion::createDispersion("%s")',data.dashD);
    end

    alias_str = [alias_str sprintf('if (boost::to_upper(name) == "%s") {\n',data.name)];
    for k = 1:length(data.weights)
        alias_str = [alias_str sprintf('    superfun->addFunctional(Functional::createFunctional("%s",npoints,deriv),%24.16E);\n',...
            data.functionals{k},data.weights(k))];
    end
    alias_str = [alias_str sprintf('    superfun.name_ = "%s";\n',data.name)];
    alias_str = [alias_str sprintf('    superfun.description = "%s";\n',data.description)];
    alias_str = [alias_str sprintf('    superfun.citation_ = "%s";\n',data.citation)];
    alias_str = [alias_str sprintf('    superfun.exact_exchange_ = %24.16E;\n',data.exact_exchange)];
    alias_str = [alias_str sprintf('    superfun.pt2_ = %24.16E;\n',data.pt2)];
    alias_str = [alias_str sprintf('    superfun.omega_ = %24.16E;\n',data.omega)];
    alias_str = [alias_str sprintf('    superfun.dashD_weight_ = %24.16E;\n',data.dashD_weight)];
    alias_str = [alias_str sprintf('    superfun.dashD_ = %s;\n',dashD_treatment)];
    alias_str = [alias_str sprintf('}\n')];
    
    name = fgetl(alias_ctrl);
end

fclose(alias_ctrl);

replaceInFile('SUPERFUNCTIONAL_ALIAS',alias_str,'factory.cc');
            

%=========== END OF SUPERFUNCTIONAL ADDITION ===============%

chdir('../')

