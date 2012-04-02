"""Module with utility functions used by several Python functions."""
import os
import pickle
import PsiMod
import input
from psiexceptions import *


def kwargs_lower(kwargs):
    """Function to rebuild and return *kwargs* dictionary
    with all keys made lowercase. Should be called by every
    function that could be called directly by the user.

    """
    caseless_kwargs = {}
    for key, value in kwargs.iteritems():
        caseless_kwargs[key.lower()] = value
    return caseless_kwargs


def get_psifile(fileno, pidspace=str(os.getpid())):
    """Function to return the full path and filename for psi file
    *fileno* (e.g., psi.32) in current namespace *pidspace*.

    """
    psioh = PsiMod.IOManager.shared_object()
    psio = PsiMod.IO.shared_object()
    filepath = psioh.get_file_path(fileno)
    namespace = psio.get_default_namespace()
    targetfile = filepath + 'psi' + '.' + pidspace + '.' + namespace + '.' + str(fileno)
    return targetfile


def format_molecule_for_input(mol):
    """Function to return a string of the output of
    :py:func:`input.process_input` applied to the XYZ
    format of molecule *mol*. Used to capture molecule
    information for distributed (sow/reap) input files.

    """
    commands = ''
    commands += 'input.process_input("""\nmolecule %s {\n' % (mol.name())
    commands += mol.save_string_xyz()
    commands += 'units angstrom\n}\n""", 0)\n'
    return eval(commands)


def format_options_for_input():
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Does not cover local (as opposed to global) options.

       - Does not work with array-type options.

    """
    commands = ''
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
    """Function to pickle to file *filename* the options dictionary
    *kwargs*. Mode *lmode* =2 pickles appropriate settings for
    reap mode. Used to capture Python options information for
    distributed (sow/reap) input files.

    """
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


def drop_duplicates(seq):
    """Function that given an array *seq*, returns an array without any duplicate
    entries. There is no guarantee of which duplicate entry is dropped.

    """
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes
