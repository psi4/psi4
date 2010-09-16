function createFunctionalLibrary(control,alias,templates,target)

ctrl = fopen(control);
alias_ctrl = fopen(alias);

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
printer_str_xc = '';
name_str = '';

disp(sprintf('Creating Functional Library in directory %s',target));


name = fgetl(ctrl);
while ischar(name)
    disp(sprintf('  Creating Functional %s',name));
    getF = inline(sprintf('get%sFunctional()',name));
    data = getF(0);
    
    createFunctional(data);
    header_str = [header_str sprintf('#include "%s_functional.h"\n',name)];
    constructor_str = [constructor_str sprintf('if (boost::to_upper_copy(name) == "%s")\n',name)];
    constructor_str = [constructor_str sprintf('    return boost::shared_ptr<Functional> (new %s_Functional(npoints,deriv));\n',name)];
    name_str = [name_str sprintf('names.push_back("%s");\n',name)];
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
    if (strcmpi(data.type, 'x'))
        printer_str_x = [printer_str_x sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    elseif (strcmpi(data.type, 'c'))
        printer_str_c = [printer_str_c sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    elseif (strcmpi(data.type, 'xc'))
        printer_str_xc = [printer_str_xc sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    end
    name = fgetl(ctrl);
end

fclose(ctrl);

replaceInFile('HEADERS',header_str,'factory.cc');
replaceInFile('FUNCTIONAL_CONSTRUCTORS',constructor_str,'factory.cc');
replaceInFile('FUNCTIONAL_NAMES',name_str,'factory.cc');
replaceInFile('X_FUNCTIONALS',printer_str_x,'factory.cc');
replaceInFile('C_FUNCTIONALS',printer_str_c,'factory.cc');
replaceInFile('XC_FUNCTIONALS',printer_str_xc,'factory.cc');

%============= END OF FUNCTIONAL ADDITION ===============%

%=========== BEGINNING OF SUPERFUNCTIONAL ADDITION ===============%

alias_str = '';
name_str = '';
printer_str_x = '';
printer_str_c = '';
printer_str_xc = '';

name = fgetl(alias_ctrl);
while ischar(name)
    disp(sprintf('  Creating SuperFunctional %s',name));
    getF = inline(sprintf('get%sSuperFunctional()',name));
    data = getF(0);

    dashD_treatment = sprintf('Dispersion::createDispersion("%s")',data.dashD);

    name_str = [name_str sprintf('names.push_back("%s");\n',data.name)];
    alias_str = [alias_str sprintf('if (boost::to_upper_copy(name) == "%s") {\n',data.name)];
    for k = 1:length(data.weights)
        alias_str = [alias_str sprintf('    superfun->addFunctional(Functional::createFunctional("%s",npoints,deriv),%24.16E);\n',...
            data.functionals{k},data.weights(k))];
    end
    alias_str = [alias_str sprintf('    superfun->setName("%s");\n',data.name)];
    alias_str = [alias_str sprintf('    superfun->setDescription("%s");\n',data.description)];
    alias_str = [alias_str sprintf('    superfun->setCitation("%s");\n',data.citation)];
    alias_str = [alias_str sprintf('    superfun->setExactExchange(%24.16E);\n',data.exact_exchange)];
    alias_str = [alias_str sprintf('    superfun->setPT2(%24.16E);\n',data.pt2)];
    alias_str = [alias_str sprintf('    superfun->setOmega(%24.16E);\n',data.omega)];
    alias_str = [alias_str sprintf('    superfun->setDashD(%s,%24.16E);\n',dashD_treatment,data.dashD_weight)];
    alias_str = [alias_str sprintf('}\n')];
    
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
    if (data.exact_exchange ~= 0.0)
        hybrid = 'X';
    else
        hybrid = ' ';
    end
    if (data.pt2 ~= 0.0)
        pt2 = 'X';
    else
        pt2 = ' ';
    end
    if (data.omega ~= 0.0)
        rc = 'X';
    else
        rc = ' ';
    end
    if (data.dashD_weight ~= 0.0)
        dashD = 'X';
    else
        dashD = ' ';
    end
    if (strcmpi(data.type, 'x'))
        printer_str_x = [printer_str_x sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
            name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
    elseif (strcmpi(data.type, 'c'))
        printer_str_c = [printer_str_c sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
            name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
    elseif (strcmpi(data.type,'xc'))
        printer_str_xc = [printer_str_xc sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
            name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
    end
    name = fgetl(alias_ctrl);
end

fclose(alias_ctrl);

replaceInFile('SUPER_CONSTRUCTORS',alias_str,'factory.cc');
replaceInFile('SUPER_NAMES',name_str,'factory.cc');
replaceInFile('X_SUPERFUNCTIONALS',printer_str_x,'factory.cc');
replaceInFile('C_SUPERFUNCTIONALS',printer_str_c,'factory.cc');
replaceInFile('XC_SUPERFUNCTIONALS',printer_str_xc,'factory.cc');
            

%=========== END OF SUPERFUNCTIONAL ADDITION ===============%

chdir('../')

