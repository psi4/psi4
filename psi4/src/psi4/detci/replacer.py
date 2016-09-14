import glob


files = glob.glob('*.cc')
files += glob.glob('*.h')

#rep_list = ['ci_tol.h', 'civect.h', 'ciwave.h', 'globaldefs.h odometer.h', 'slaterd.h', 'structs.h']
rep_list = ['globaldefs.h', 'odometer.h']

for fn in files:
    with open(fn, 'r') as f:
        data = f.read()

    for string in rep_list:    
        fname = '#include "**REP**"'
        rname = '#include "psi4/detci/**REP**"'
        fname = fname.replace('**REP**', string)
        rname = rname.replace('**REP**', string)
        data = data.replace(fname, rname)

    with open(fn, 'w') as f:
        f.write(data)


