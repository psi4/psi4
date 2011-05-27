function [lib] = getDefaultLibraryData()

% Absolute path to library template files
lib.templates = '/theoryfs2/ds/parrish/matlab/templates/';
% Absolute path to object directory of individual .h / .cc files
lib.objdir = '/theoryfs2/ds/parrish/matlab/objdir2/';
% Absolute path of installation of final library
lib.installdir = '/theoryfs2/ds/parrish/matlab/install/';
% Absolute path for test file outputs
lib.testdir = '/theoryfs2/ds/parrish/matlab/tests/';
% Test unmade/remake functionals in matlab?
lib.test = false;

lib.compile = true;
lib.install = true;

% List of all functionals
lib.functionals = {...
    'S_X',
    'B_X',
    'B88_X',
    'B88II_X',
    'PBE_X',
    'PW91_X',
    'FT97B_X',
    'LYP_C',
    'VWN5RPA_C',
    'VWN5_C',
    'PZ81_C',
    'P86_C'
    'PW91_C',
    'PW92_C',
    'PBE_C',
    'FT97_C',
    'EDF1',
    'EDF1_X',
    'EDF1_C',
    'B97_0',
    'B97_1',
    'B97_2',
    'B97_D2',
    'HCTH',
    'HCTH407',
    'HCTH147',
    'HCTH120',
    'M05',
    'M05_2X',
    'TauHCTH',
    'TauHCTH0',
    'wB97',
    'wB97X'
};

% List of all functionals to remake
lib.remakes = {
};

% Superfunctional file
lib.superfunctionals = '/theoryfs2/ds/parrish/matlab/superfunctionals.xls'; 
lib.superfunctionals_sheet = 'SuperFunctionals'; 


