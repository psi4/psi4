function buildFunctional(name, install, paths)

if (nargin < 2)
    install = true;
end

if (nargin < 3)
    paths.definition = '.';
    paths.install = '../../';
    paths.python = '../mfiles/';
    paths.templates = '../templates/';
    paths.driver = '../../../../../lib/python/';
end

disp(sprintf('\n'))
disp('      *****************************************')
disp('      *                                       *')
disp('      *    Autofunctional Generator v. 2.0    *')
disp('      *              Rob Parrish              *')
disp('      *          The Sherrill Group           *')
disp('      *    Georgia Institute of Technology    *')
disp('      *                                       *')
disp('      *****************************************')
disp(sprintf('\n'))
disp('')

disp(sprintf('Creating Functional      %s', name))
disp(sprintf('Compiling Functional in  %s', paths.definition))
if (install)
    disp(sprintf('Installing Functional in %s', paths.install))
end

paths.current = pwd();
cd(paths.definition)
system(sprintf('/bin/cp %s/*py .', paths.python));

getF = inline(sprintf('get%sFunctional()',name));
functional = getF(0);

% => Preamble <= %

fh = fopen('preamble', 'w');
fprintf(fh, 'name_ = "%s";\n', functional.name);
fprintf(fh, 'description_ = "    %s\\n";\n', functional.description);
fprintf(fh, 'citation_ = "    %s\\n";\n', functional.citation);
fprintf(fh, 'alpha_ = 1.0;\n');
fprintf(fh, 'omega_ = 0.0;\n');
fprintf(fh, 'lrc_ = false;\n');
if (functional.is_gga)
    fprintf(fh, 'gga_ = true;\n');
else
    fprintf(fh, 'gga_ = false;\n');
end    
if (functional.is_meta)
    fprintf(fh, 'meta_ = true;\n');
else
    fprintf(fh, 'meta_ = false;\n');
end    
for k = 1:length(functional.param_names)
    fprintf(fh, 'parameters_["%s"] = %24.16E;\n', functional.param_names{k}, functional.param_vals(k));
end

fclose(fh);

% => Parameters <= %

fh = fopen('parameters', 'w');
for k = 1:length(functional.param_names)
    fprintf(fh, 'double %s = parameters_["%s"];\n', functional.param_names{k}, functional.param_names{k});
end
fclose(fh);

% => Partials <= %

deriv = 2;
if (functional.deriv2 == 0)
    deriv = 1;
end

buildPartials(functional.functional, 'functional', deriv)
buildPartials(functional.functional_a0, 'functional_rho_a0', deriv)
buildPartials(functional.functional_b0, 'functional_rho_b0', deriv)

system(sprintf('/bin/cp %s/functional.cc %sfunctional.cc', paths.templates, functional.name));
system(sprintf('/bin/cp %s/functional.h %sfunctional.h', paths.templates, functional.name));
system(sprintf('./assemble.py %s', functional.name));

if (install)
    system(sprintf('./install.py %s %s %s', functional.name, paths.install, paths.driver));
end

system(sprintf('rm *py'));

cd(paths.current)
