function [lib] = getDefaultLibraryDataRev()

% Absolute path to library template files
lib.templates = '/theoryfs2/ds/parrish/matlab/templates_rev/';
% Absolute path to object directory of individual .h / .cc files
lib.objdir = '/theoryfs2/ds/parrish/matlab/objdir_rev/';
% Absolute path of installation of final library
lib.installdir = '/theoryfs2/ds/parrish/matlab/install_rev/';
% Absolute path for test file outputs
lib.testdir = '/theoryfs2/ds/parrish/matlab/tests_rev/';
% Test unmade/remake functionals in matlab?
lib.test = false;

lib.compile = true;
lib.install = true;

% List of all functionals
lib.functionals = {...
    'SII_X',
    'B88II_X',
    'BII_X',
    'BX_X',
    'PW92II_C',
    'PBEII_C'
};

% List of all functionals to remake
lib.remakes = {
};

% Superfunctional file
lib.superfunctionals = '/theoryfs2/ds/parrish/matlab/superfunctionals_rev.xls'; 
lib.superfunctionals_sheet = 'SuperFunctionals'; 


