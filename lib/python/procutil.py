import os
import pickle
import PsiMod
import input
from psiexceptions import *


# Re-builds kwargs with only lowercase keyword names. Should be called by every
#   function that could be called directly by the user
def kwargs_lower(kwargs):
    caseless_kwargs = {}
    for key, value in kwargs.iteritems():
        caseless_kwargs[key.lower()] = value
    return caseless_kwargs

def get_psifile(fileno, pidspace=str(os.getpid())):
    psioh = PsiMod.IOManager.shared_object()
    psio = PsiMod.IO.shared_object()
    filepath = psioh.get_file_path(fileno)
    namespace = psio.get_default_namespace()
    targetfile =  filepath + 'psi' + '.' + pidspace + '.' + namespace + '.' + str(fileno)
    return targetfile

# Common sow/reap functions
def format_molecule_for_input(mol):
    commands  = ''
    commands += 'input.process_input("""\nmolecule %s {\n' % (mol.name())
    commands += mol.save_string_xyz()
    commands += 'units angstrom\n}\n""")\n'
    return eval(commands)

def format_options_for_input():
    # build string of commands for options from the input file  TODO: handle local options too
    commands  = ''
    commands += """\nPsiMod.set_memory(%s)\n\n""" % (PsiMod.get_memory())
    for chgdopt in PsiMod.get_global_option_list():
        if PsiMod.has_option_changed(chgdopt):
            chgdoptval = PsiMod.get_global_option(chgdopt)
            if isinstance(chgdoptval, basestring):
                commands += """PsiMod.set_global_option('%s', '%s')\n""" % (chgdopt, chgdoptval)
            elif isinstance(chgdoptval, int) or isinstance(chgdoptval, float):
                commands += """PsiMod.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
            else:
                raise ValidationError('Option \'%s\' is not of a type (string, int, float, bool) that can be processed.' % (chgdopt))
    return commands

def format_kwargs_for_input(filename, lmode=1, **kwargs):
    if lmode == 2:
        kwargs['mode'] = 'reap'
        kwargs['linkage'] = os.getpid()
    filename.write('''\npickle_kw = ("""''')
    pickle.dump(kwargs, filename)
    filename.write('''""")\n''')
    filename.write("""\nkwargs = pickle.loads(pickle_kw)\n""")
    if lmode == 2:
        kwargs['mode'] = 'sow'
        del kwargs['linkage']

