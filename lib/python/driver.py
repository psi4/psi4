import PsiMod
from proc import *

#Procedure lookup tables
procedures = {'energy' : {
            'scf'           : run_scf,
            'hf'            : run_scf,
            'rhf'           : run_scf,
            'uhf'           : run_scf,
            'rohf'          : run_scf,
            'rks'           : run_scf,
            'uks'           : run_scf,
            'dfmp2'         : run_dfmp2,
            'scs-dfmp2'     : run_dfmp2,
            'sapt0'         : run_sapt,
            'scs-sapt'      : run_sapt,
            'sapt2'         : run_sapt,
            'sapt2+'        : run_sapt,
            'sapt2+3'       : run_sapt,
            'sapt_dft'      : run_sapt,
            'ccsd'          : run_ccsd,
            'ccsd(t)'       : run_ccsd_t
        },
        'gradient' : {
        },
        'hessian' : { }}

def energy(name, **kwargs):

    if (kwargs.has_key('molecule')):
        activate(kwargs['molecule'])
    if (kwargs.has_key('bases')):
        pass
    if (kwargs.has_key('functional')):
        pass

    try:
        return procedures['energy'][name](name,**kwargs)
    except KeyError:
        raise SystemExit('Energy Method %s Not Defined' %(name))

def gradient(name, **kwargs):
    try:
        return procedures['gradient'][name](name, **kwargs)
    except KeyError:
        raise SystemExit('Gradient Method %s Not Defined' %(name))

def hessian(name, **kwargs):
    try:
        return procedures['hessian'][name](name, **kwargs)
    except KeyError:
        raise SystemExit('Hessian Method %s Not Defined' %(name))


