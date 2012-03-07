function createFunctionalLibrary(libdata)

%=========== STARTUP ===============%

% Grab some handles
if (nargin < 1)
    libdata = getDefaultLibraryData();
end

original_dir = pwd;
templates = libdata.templates;
objdir = libdata.objdir;
installdir = libdata.installdir;
testdir = libdata.testdir;
test = libdata.test;
superfunctionals = libdata.superfunctionals;

% Titles and paths
disp(sprintf('\n'))
disp('      *****************************************')
disp('      *                                       *')
disp('      *    Autofunctional Generator v. 1.0    *')
disp('      *              Rob Parrish              *')
disp('      *          The Sherrill Group           *')
disp('      *    Georgia Institute of Technology    *')
disp('      *                                       *')
disp('      *****************************************')
disp(sprintf('\n'))
disp('')

disp(sprintf('Taking templates from:     %s',templates));
disp(sprintf('Taking aliases from:       %s',superfunctionals));
disp(sprintf('Compiling in:              %s',objdir));
disp(sprintf('Installing to:             %s',installdir));
if (test)
    disp(sprintf('Testing in:                %s',testdir));
end
disp(sprintf('\n'))

%=========== END STARTUP ===============%

%=========== BEGINNING OF FUNCTIONAL COMPILATION ===============%

% total strings for factory.cc
header_str = '';
constructor_str = '';
printer_str_x = '';
printer_str_c = '';
printer_str_xc = '';
name_str = '';

if (~isdir(objdir))
    mkdir(objdir);
end

for k = 1:length(libdata.functionals)
    name = libdata.functionals{k};
    cd(objdir);
    
    % Does the dir exist and carry the file metadata.m?
    % And the .cc and .h files?
    dir_exists = isdir(name);
    metadata_exists = false;
    cpp_exists = false;
    if (dir_exists)
       this_obj = dir(name);
       metadata_exists = any(strcmpi('metadata.mat',{this_obj.name})); 
       cpp_exists = any(strcmpi([name '_functional.cc'],{this_obj.name})) && ...
            any(strcmpi([name '_functional.h'],{this_obj.name})); 
    end    

    if (~any(strcmpi(name,libdata.remakes)) && metadata_exists && cpp_exists)
        % Metadata exists and user does not wish a remake
        % Grab the metadata
        cd(name);
        load metadata this_header_str this_constructor_str this_printer_str_x this_printer_str_c this_printer_str_xc this_name_str;
        header_str = [header_str this_header_str];
        constructor_str = [constructor_str this_constructor_str];
        printer_str_x = [printer_str_x this_printer_str_x];
        printer_str_c = [printer_str_c this_printer_str_c];
        printer_str_xc = [printer_str_xc this_printer_str_xc];
        name_str = [name_str this_name_str];
        disp(sprintf('  Creating Functional %s ... premade.',name));
        continue;
    end

    % If you're here, you need to make the functional
    if (any(strcmpi(name,libdata.remakes))) 
        disp(sprintf('  Creating Functional %s ... remake.',name));
    else
        disp(sprintf('  Creating Functional %s ... new make.',name));
    end   
 
    % Grab the functional
    getF = inline(sprintf('get%sFunctional()',name));
    data = getF(0);

    % Test functional if needed
    if (test)
        chdir(testdir)
        testFunctional(data,1,[data.name '_Test.dat']);
        chdir(objdir)
    end

    % Create the compilation folder   
    if (dir_exists)
        system(sprintf('rm -rf %s',name));
    end
    mkdir(name);

    % Go into the compilation folder: 
    chdir(name);    

    % Put the templates into the compilation folder
    system(sprintf('cp %sfunctional.cc.template functional.cc.template',templates));
    system(sprintf('cp %sfunctional.h.template functional.h.template',templates));

    % Call the build subroutine
    if (data.algorithm == 1)
        createFunctional(data);
    elseif (data.algorithm == 2)
        createFunctionalII(data);
    end   

    % Get the metadata 
    this_header_str = [sprintf('#include "%s_functional.h"\n',name)];
    this_constructor_str = [sprintf('if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("%s")))\n',name)];
    this_constructor_str = [this_constructor_str sprintf('    return boost::shared_ptr<Functional> (new %s_Functional(npoints,deriv));\n',name)];
    this_name_str = [sprintf('names.push_back("%s");\n',name)];
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
    this_printer_str_x = '';
    this_printer_str_c = '';
    this_printer_str_xc = '';
    if (strcmpi(data.type, 'x'))
        this_printer_str_x = [sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    elseif (strcmpi(data.type, 'c'))
        this_printer_str_c = [sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    elseif (strcmpi(data.type, 'xc'))
        this_printer_str_xc = [sprintf('f << "   %-10s      %s         %s         %s " << endl;\n',name,lsda,gga,meta)];
    end
    
    %Add the metadata to the total metadata strings
    header_str = [header_str this_header_str];
    constructor_str = [constructor_str this_constructor_str];
    printer_str_x = [printer_str_x this_printer_str_x];
    printer_str_c = [printer_str_c this_printer_str_c];
    printer_str_xc = [printer_str_xc this_printer_str_xc];
    name_str = [name_str this_name_str];

    %Save the metadata 
    save metadata this_header_str this_constructor_str this_printer_str_x this_printer_str_c this_printer_str_xc this_name_str;
    
end

%============= END OF FUNCTIONAL COMPILATION ===============%

%=========== BEGINNING OF SUPERFUNCTIONAL COMPILATION ===============%
if (libdata.install)

    % All superfunctionals are remade on every pass, this does not
    % require much work at all
   
    % Grab the superfunctionals
    supers = getSuperFunctionals(superfunctionals);
 
    % Total strings for factory.cc
    super_alias_str = '';
    super_name_str = '';
    super_printer_str_x = '';
    super_printer_str_c = '';
    super_printer_str_xc = '';
    
    for k = 1:length(supers)
    
        data = supers(k);
        name = data.name;
        disp(sprintf('  Creating SuperFunctional %s',name));
    
        
        if (strcmpi(data.dashD, 'none')) 
            dashD_treatment = '';
        elseif (strcmpi(data.dashD, '-D1'))
            dashD_treatment = sprintf('    superfun->setDashD(Dispersion::createDispersion("%s",%24.16E),1.0);\n','-D1',data.dashD_s6);
        elseif (strcmpi(data.dashD, '-D2'))
            dashD_treatment = sprintf('    superfun->setDashD(Dispersion::createDispersion("%s",%24.16E),1.0);\n','-D2',data.dashD_s6);
        end   
 
        super_name_str = [super_name_str sprintf('names.push_back("%s");\n',data.name)];
        super_alias_str = [super_alias_str sprintf('if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("%s"))) {\n',data.name)];
        for k = 1:length(data.weights)
            super_alias_str = [super_alias_str sprintf('    superfun->addFunctional(Functional::createFunctional("%s",npoints,deriv),%24.16E);\n',...
                data.functionals{k},data.weights(k))];
        end
        super_alias_str = [super_alias_str sprintf('    superfun->setName("%s");\n',data.name)];
        super_alias_str = [super_alias_str sprintf('    superfun->setDescription("%s");\n',data.description)];
        super_alias_str = [super_alias_str sprintf('    superfun->setCitation("%s");\n',data.citation)];
        super_alias_str = [super_alias_str sprintf('    superfun->setExactExchange(%24.16E);\n',data.exact_exchange)];
        super_alias_str = [super_alias_str sprintf('    superfun->setPT2(%24.16E);\n',data.pt2)];
        super_alias_str = [super_alias_str sprintf('    superfun->setOmega(%24.16E);\n',data.omega)];
        % -D treatment
        super_alias_str = [super_alias_str dashD_treatment];
    
        super_alias_str = [super_alias_str sprintf('}\n')];
        
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
        if (data.dashD_s6 ~= 0.0)
            dashD = 'X';
        else
            dashD = ' ';
        end
        if (strcmpi(data.type, 'x'))
            super_printer_str_x = [super_printer_str_x sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
                name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
        elseif (strcmpi(data.type, 'c'))
            super_printer_str_c = [super_printer_str_c sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
                name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
        elseif (strcmpi(data.type,'xc'))
            super_printer_str_xc = [super_printer_str_xc sprintf('f << "   %-10s      %s         %s         %s         %s         %s         %s         %s " << endl;\n',...
                name,lsda,gga,meta,hybrid,pt2,rc,dashD)];
        end
    end
    
end
%=========== END OF SUPERFUNCTIONAL COMPILATION ===============%

%=========== INSTALLATION ===============%
if (libdata.install)
    
    % Wipe the install dir and remake
    if (isdir(installdir))
        system(sprintf('rm -rf %s', installdir));
    end
    mkdir(installdir);
    
    % Place functionals in the install dir
    for k = 1:length(libdata.functionals)
        name = libdata.functionals{k};
        system(sprintf('cp %s%s/*.cc %s',objdir, name, installdir));
        system(sprintf('cp %s%s/*.h %s',objdir, name, installdir));
    end
    
    % Move all cc, h, and in files from templates to install
    system(sprintf('cp %s*.cc %s', templates,installdir));
    system(sprintf('cp %s*.h %s', templates,installdir));
    system(sprintf('cp %s*.in %s', templates,installdir));
    
    % Copy factory.cc.template to desired $INSTALL/factory.cc file
    install_dir = dir(installdir);
    if (any(strcmpi(['factory.cc'],{install_dir.name})))
        system(sprintf('rm "factory.cc"'))
    end
    system(sprintf('cp %sfactory.cc.template %sfactory.cc',templates, installdir));
    factory_path = [installdir 'factory.cc'];
    
    % factory.cc for functionals
    replaceInFile('HEADERS',header_str,factory_path);
    replaceInFile('FUNCTIONAL_CONSTRUCTORS',constructor_str,factory_path);
    replaceInFile('FUNCTIONAL_NAMES',name_str,factory_path);
    replaceInFile('X_FUNCTIONALS',printer_str_x,factory_path);
    replaceInFile('XC_FUNCTIONALS',printer_str_xc,factory_path);
    replaceInFile('C_FUNCTIONALS',printer_str_c,factory_path);
    
    % factory.cc for superfunctionals
    replaceInFile('SUPER_CONSTRUCTORS',super_alias_str,factory_path);
    replaceInFile('SUPER_NAMES',super_name_str,factory_path);
    replaceInFile('X_SUPERFUNCTIONALS',super_printer_str_x,factory_path);
    replaceInFile('XC_SUPERFUNCTIONALS',super_printer_str_xc,factory_path);
    replaceInFile('C_SUPERFUNCTIONALS',super_printer_str_c,factory_path);

    % Copy templates to install
    mkdir(sprintf('%s/templates/',installdir));
    system(sprintf('cp %s/* %s/templates/', templates,installdir));
    % Copy matlab to install
    mkdir(sprintf('%s/matlab/',installdir));
    system(sprintf('cp %s/*.m %s/matlab/', original_dir, installdir));
end            
%=========== END INSTALLATION ===============%

chdir(original_dir)

%=========== END ===============%

end
