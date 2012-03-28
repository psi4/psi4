function buildLibrary(list)

if (nargin < 1)
    list = { 'FT97B_X', 
        'PZ81_C',
        'P86_C',
        'VWN_C',
        'LYP_C',
        'PW91_C',
        'PW92_C',
        'PBE_C',
        'FT97_C',
        'B972_C',
        'B974_C'
    };
end

current = pwd();

for k = 1:length(list)
    name = list{k};
    cd(name)
    buildFunctional(name)
    cd('../')
end 

