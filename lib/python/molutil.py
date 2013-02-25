"""Module with utility functions that act on molecule objects."""
import os
import re
import subprocess
import socket
import shutil
import random
import math
import PsiMod
import physconst
from text import *
from dashparam import *


def extract_clusters(mol, ghost=True, cluster_size=0):
    """Function to return all subclusters of the molecule *mol* of
    real size *cluster_size* and all other atoms ghosted if *ghost*
    equals true, all other atoms discarded if *ghost* is false. If
    *cluster_size* = 0, returns all possible combinations of cluster size.

    """
    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []
    ghosts = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            if (ghost):
                ghosts = []
                for g in range(nfrag, 0, -1):
                    if (g not in reals):
                        ghosts.append(g)
                #print "Cluster #%d: %s reals, %s ghosts" % (counter,str(reals), str(ghosts))
                clusters.append(mol.extract_subsets(reals, ghosts))
            else:
                #print "Cluster #%d: %s reals" % (counter,str(reals))
                clusters.append(mol.extract_subsets(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


def extract_cluster_indexing(mol, cluster_size=0):
    """Function to returns a LIST of all subclusters of the molecule *mol* of
    real size *cluster_size*. If *cluster_size* = 0, returns all possible
    combinations of cluster size.

    """
    import copy

    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            clusters.append(copy.deepcopy(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


def new_set_attr(self, name, value):
    """Function to redefine set_attr method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)
    if isvar:
        fxn = object.__getattribute__(self, "set_variable")
        fxn(name, value)
        return

    object.__setattr__(self, name, value)


def new_get_attr(self, name):
    """Function to redefine get_attr method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)

    if isvar:
        fxn = object.__getattribute__(self, "get_variable")
        return fxn(name)

    return object.__getattribute__(self, name)


def BFS(self):
    """Perform a breadth-first search (BFS) on the real atoms
    in molecule, returning an array of atom indices of fragments.
    Relies upon van der Waals radii and so faulty for close
    (esp. hydrogen-bonded) fragments. Original code from
    Michael S. Marshall.

    """
    vdW_diameter = {
        'H':  1.001 / 1.5,
        'HE': 1.012 / 1.5,
        'LI': 0.825 / 1.5,
        'BE': 1.408 / 1.5,
        'B':  1.485 / 1.5,
        'C':  1.452 / 1.5,
        'N':  1.397 / 1.5,
        'O':  1.342 / 1.5,
        'F':  1.287 / 1.5,
        'NE': 1.243 / 1.5,
        'NA': 1.144 / 1.5,
        'MG': 1.364 / 1.5,
        'AL': 1.639 / 1.5,
        'SI': 1.716 / 1.5,
        'P':  1.705 / 1.5,
        'S':  1.683 / 1.5,
        'CL': 1.639 / 1.5,
        'AR': 1.595 / 1.5}

    Queue = []
    White = range(self.natom())  # untouched
    Black = []  # touched and all edges discovered
    Fragment = []  # stores fragments

    start = 0  # starts with the first atom in the list
    Queue.append(start)
    White.remove(start)

    # Simply start with the first atom, do a BFS when done, go to any
    #   untouched atom and start again iterate until all atoms belong
    #   to a fragment group
    while len(White) > 0 or len(Queue) > 0:  # Iterates to the next fragment
        Fragment.append([])

        while len(Queue) > 0:                # BFS within a fragment
            for u in Queue:                  # find all (still white) nearest neighbors to vertex u
                for i in White:
                    dist = physconst.psi_bohr2angstroms * math.sqrt((self.x(i) - self.x(u)) ** 2 + \
                        (self.y(i) - self.y(u)) ** 2 + (self.z(i) - self.z(u)) ** 2)
                    if dist < vdW_diameter[self.symbol(u)] + vdW_diameter[self.symbol(i)]:
                        Queue.append(i)      # if you find you, put in the queue
                        White.remove(i)      # and remove it from the untouched list
            Queue.remove(u)                  # remove focus from Queue
            Black.append(u)
            Fragment[-1].append(int(u))      # add to group (0-indexed)
            Fragment[-1].sort()              # preserve original atom ordering

        if len(White) != 0:                  # can't move White -> Queue if no more exist
            Queue.append(White[0])
            White.remove(White[0])

    return Fragment


def run_dftd3(self, func=None, dashlvl=None, dashparam=None, dertype=None):
    """Function to call Grimme's dftd3 program (http://toc.uni-muenster.de/DFTD3/)
    to compute the -D correction of level *dashlvl* using parameters for
    the functional *func*. The dictionary *dashparam* can be used to supply
    a full set of dispersion parameters in the absense of *func* or to supply
    individual overrides in the presence of *func*. Returns energy if *dertype* is 0,
    gradient if *dertype* is 1, else tuple of energy and gradient if *dertype*
    unspecified. The dftd3 executable must be independently compiled and found in 
    :envvar:`PATH`.

    """
    # Validate arguments
    if self is None:
        self = PsiMod.get_active_molecule()

    dashlvl = dashlvl.lower()
    dashlvl = dash_alias['-' + dashlvl][1:] if ('-' + dashlvl) in dash_alias.keys() else dashlvl
    if dashlvl not in dashcoeff.keys():
        raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))

    if dertype is None:
        dertype = -1
    elif der0th.match(str(dertype)):
        dertype = 0
    elif der1st.match(str(dertype)):
        dertype = 1
    elif der2nd.match(str(dertype)):
        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))
    else:
        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))

    if func is None:
        if dashparam is None:
            # defunct case
            raise ValidationError("""Parameters for -D correction missing. Provide a func or a dashparam kwarg.""")
        else:
            # case where all param read from dashparam dict (which must have all correct keys)
            func = 'custom'
            dashcoeff[dashlvl][func] = {}
            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
            for key in dashcoeff[dashlvl]['b3lyp'].keys():
                if key in dashparam.keys():
                    dashcoeff[dashlvl][func][key] = dashparam[key]
                else:
                    raise ValidationError("""Parameter %s is missing from dashparam dict %s.""" % (key, dashparam))
    else:
        func = func.lower()
        if func not in dashcoeff[dashlvl].keys():
            raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))
        if dashparam is None:
            # (normal) case where all param taken from dashcoeff above
            pass
        else:
            # case where items in dashparam dict can override param taken from dashcoeff above
            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
            for key in dashcoeff[dashlvl]['b3lyp'].keys():
                if key in dashparam.keys():
                    dashcoeff[dashlvl][func][key] = dashparam[key]

    # Move ~/.dftd3par.<hostname> out of the way so it won't interfere
    defaultfile = os.path.expanduser('~') + '/.dftd3par.' + socket.gethostname()
    defmoved = False
    if os.path.isfile(defaultfile):
        os.rename(defaultfile, defaultfile + '_hide')
        defmoved = True

    # Setup unique scratch directory and move in
    current_directory = os.getcwd()
    psioh = PsiMod.IOManager.shared_object()
    psio = PsiMod.IO.shared_object()
    os.chdir(psioh.get_default_path())
    dftd3_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
        '.dftd3.' + str(random.randint(0, 99999))
    if os.path.exists(dftd3_tmpdir) is False:
        os.mkdir(dftd3_tmpdir)
    os.chdir(dftd3_tmpdir)

    # Write dftd3_parameters file that governs dispersion calc
    paramfile = './dftd3_parameters'
    pfile = open(paramfile, 'w')
    pfile.write(dash_server(func, dashlvl, 'dftd3'))
    pfile.close()

    # Write dftd3_geometry file that supplies geometry to dispersion calc
    geomfile = './dftd3_geometry.xyz'
    gfile = open(geomfile, 'w')
    numAtoms = self.natom()
    geom = self.save_string_xyz()
    reals = []
    for line in geom.splitlines():
      if line.split()[0] == 'Gh':
        numAtoms -= 1
      else:
        reals.append(line)
        
    gfile.write(str(numAtoms)+'\n')
    for line in reals:
      gfile.write(line.strip()+'\n')
    gfile.close()

    # Call dftd3 program
    try:
        dashout = subprocess.Popen(['dftd3', geomfile, '-grad'], stdout=subprocess.PIPE)
    except OSError:
        raise ValidationError('Program dftd3 not found in path.')
    out, err = dashout.communicate()

    # Parse output (could go further and break into E6, E8, E10 and Cn coeff)
    success = False
    for line in out.splitlines():
        if re.match(' Edisp /kcal,au', line):
            sline = line.split()
            dashd = float(sline[3])
        if re.match(' normal termination of dftd3', line):
            success = True

    if not success:
        raise ValidationError('Program dftd3 did not complete successfully.')

    # Parse grad output
    derivfile = './dftd3_gradient'
    dfile = open(derivfile, 'r')
    dashdderiv = []
    i = 0
    for line in geom.splitlines():
      if i == 0:
        i += 1
      else:
        if line.split()[0] == 'Gh':
          dashdderiv.append([0.0, 0.0, 0.0])
        else:
          temp = dfile.readline()
          dashdderiv.append([float(x.replace('D', 'E')) for x in temp.split()])
    dfile.close()

    if len(dashdderiv) != self.natom():
        raise ValidationError('Program dftd3 gradient file has %d atoms- %d expected.' % \
            (len(dashdderiv), self.natom()))
    psi_dashdderiv = PsiMod.Matrix(self.natom(), 3)
    psi_dashdderiv.set(dashdderiv)

    # Print program output to file if verbose
    verbose = PsiMod.get_option('SCF', 'PRINT')
    if verbose >= 3:
        PsiMod.print_out('\n  ==> DFTD3 Output <==\n')
        PsiMod.print_out(out)
        dfile = open(derivfile, 'r')
        PsiMod.print_out(dfile.read().replace('D', 'E'))
        dfile.close()
        PsiMod.print_out('\n')

    # Clean up files and remove scratch directory
    os.unlink(paramfile)
    os.unlink(geomfile)
    os.unlink(derivfile)
    if defmoved is True:
        os.rename(defaultfile + '_hide', defaultfile)

    os.chdir('..')
    try:
        shutil.rmtree(dftd3_tmpdir)
    except OSError as e:
        ValidationError('Unable to remove dftd3 temporary directory %s' % e, file=sys.stderr)
    os.chdir(current_directory)

    # return -D & d(-D)/dx
    PsiMod.set_variable('DISPERSION CORRECTION ENERGY', dashd)
    if dertype == -1:
        return dashd, dashdderiv
    elif dertype == 0:
        return dashd
    elif dertype == 1:
        return psi_dashdderiv


def dynamic_variable_bind(cls):
    """Function to dynamically add extra members to
    the PsiMod.Molecule class.

    """
    cls.__setattr__ = new_set_attr
    cls.__getattr__ = new_get_attr
    cls.BFS = BFS
    cls.run_dftd3 = run_dftd3


dynamic_variable_bind(PsiMod.Molecule)  # pass class type, not class instance


#
# Define geometry to be used by PSI4.
# The molecule created by this will be set in options.
#
# geometry("
#   O  1.0 0.0 0.0
#   H  0.0 1.0 0.0
#   H  0.0 0.0 0.0
#
def geometry(geom, name="default"):
    """Function to create a molecule object of name *name*
    from the geometry in string *geom*.

    """
    molecule = PsiMod.Molecule.create_molecule_from_string(geom)
    molecule.set_name(name)

    activate(molecule)

    return molecule


def activate(mol):
    """Function to set molecule object *mol* as the current active molecule."""
    PsiMod.set_active_molecule(mol)
    #PsiMod.IO.set_default_namespace(mol.get_name())
