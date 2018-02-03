import re

import numpy as np

from ..util import distance_matrix
from ..exceptions import *
from .chgmult import validate_and_fill_chgmult
from .nucleus import reconcile_nucleus
from .regex import *

try:
    long(1)
except NameError:
    long = int


def from_arrays(geom=None,
                elea=None,
                elez=None,
                elem=None,
                mass=None,
                real=None,
                elbl=None,
                name=None,
                units='Angstrom',
                input_units_to_au=None,
                fix_com=False,
                fix_orientation=False,
                fix_symmetry=None,
                fragment_separators=None,
                #fragment_types=None,  # keep?
                fragment_charges=None,
                fragment_multiplicities=None,
                molecular_charge=None,
                molecular_multiplicity=None,
                speclabel=True,
                tooclose=0.1,
                zero_ghost_fragments=False,
                nonphysical=False,
                mtol=1.e-3,
                verbose=1):
    """Compose a Molecule dict from unvalidated arrays and variables, returning dict.

    minimum is geom and one of elem, elez, elbl


    Parameters
    ----------
    See fields of return  molrec below. Required parameters are `geom` and one of `elem`, `elez`, `elbl` (`speclabel=True`)
    geom : array-like
        (nat, 3) or (3 * nat, ) ndarray or list o'lists of Cartesian coordinates.
    fragment_separators : array-like of int, optional
        (nfr - 1, ) list of atom indices at which to split `geom` into fragments.
    #fragment_types : array-like of {'Real', 'Ghost', 'Absent'}, optional
    #    (nfr, ) list of fragment types. If given, size must be consistent with `fragment_separators`.
    elbl : ndarray of str
        (nat, ) Label extending `elem` symbol, possibly conveying ghosting, isotope, mass, tagging information.

    tooclose : float, optional
        Interatom distance (native `geom` units) nearer than which atoms not allowed.    
    nonphysical : bool, optional

    speclabel : bool, optional
        If `True`, interpret `elbl` as potentially full nucleus spec including
        ghosting, isotope, mass, tagging information, e.g., `@13C_mine` or
        `He4@4.01`. If `False`, interpret `elbl` as only the user/tagging
        extension to nucleus label, e.g. `_mine` or `4` in the previous examples.


    Returns
    -------
    molrec : dict
        Molecule dictionary spec follows. Its principles are (1)
        contents are fully validated and defaulted - no error checking
        necessary, (2) contents may be mildly redundant - atomic
        numbers and element symbols present, (3) big system,
        nat-length single-type arrays, not small system, nat-number
        heterogeneous objects, (4) some fields are optional (e.g.,
        symmetry) but largely self-describing so units or fix_com must
        be present.

    name : str, optional
        Label for molecule; should be valid Python identifier.
    units : {'Angstrom', 'Bohr'}
        Units for `geom`.
    input_units_to_au : float, optional
        If `units='Angstrom'`, overrides consumer's value for [A]-->[a0] conversion.
    fix_com : bool
        Whether translation of `geom` is allowed or disallowed.
    fix_orientation : bool
        Whether rotation of `geom` is allowed or disallowed.
    fix_symmetry : str, optional
        Maximal point group symmetry which `geom` should be treated. Lowercase.
    geom : ndarray of float
        (3 * nat, ) Cartesian coordinates in `units`.
    elea : ndarray of int
        (nat, ) Mass number for atoms, if known isotope, else -1.
    elez : ndarray of int
        (nat, ) Number of protons, nuclear charge for atoms.
    elem : ndarray of str
        (nat, ) Element symbol for atoms.
    mass : ndarray of float
        (nat, ) Atomic mass [u] for atoms.
    real : ndarray of bool
        (nat, ) Real/ghostedness for atoms.
    elbl : ndarray of str
        (nat, ) Label with any tagging information from element spec.
    fragment_separators : list of int
        (nfr - 1, ) list of atom indices at which to split `geom` into fragments.
    fragment_charges : list of float
        (nfr, ) list of charge allocated to each fragment.
    fragment_multiplicities : list of int
        (nfr, ) list of multiplicity allocated to each fragment.
    molecular_charge : float
        total charge on system.
    molecular_multiplicity : int
        total multiplicity on system.

    """
    molinit = {}
    molinit.update(validate_and_fill_universals(name=name,
                                                units=units,
                                                input_units_to_au=input_units_to_au,
                                                fix_com=fix_com,
                                                fix_orientation=fix_orientation,
                                                fix_symmetry=fix_symmetry))

    molinit.update(validate_and_fill_geometry(  geom=geom,
                                                tooclose=tooclose))
    nat = molinit['geom'].shape[0] // 3

    molinit.update(validate_and_fill_nuclei(    nat,
                                                elea=elea,
                                                elez=elez,
                                                elem=elem,
                                                mass=mass,
                                                real=real,
                                                elbl=elbl,
                                                speclabel=speclabel,
                                                nonphysical=nonphysical,
                                                mtol=mtol,
                                                verbose=verbose))

    molinit.update(validate_and_fill_fragments( nat,
                                                fragment_separators=fragment_separators,
                                                #fragment_types=fragment_types,
                                                fragment_charges=fragment_charges,
                                                fragment_multiplicities=fragment_multiplicities))

    Z_available = molinit['elez'] * molinit['real'] * 1.
    molinit.update(validate_and_fill_chgmult(   zeff=Z_available,
                                                fragment_separators=molinit['fragment_separators'],
                                                molecular_charge=molecular_charge,
                                                fragment_charges=molinit['fragment_charges'],
                                                molecular_multiplicity=molecular_multiplicity,
                                                fragment_multiplicities=molinit['fragment_multiplicities'],
                                                zero_ghost_fragments=zero_ghost_fragments,
                                                verbose=verbose))

    return molinit


def validate_and_fill_universals(name=None,
                                 units='Angstrom',
                                 input_units_to_au=None,
                                 fix_com=False,
                                 fix_orientation=False,
                                 fix_symmetry=None):
        molinit = {}

        if name is not None:
            molinit['name'] = name

        if units in ['Angstrom', 'Bohr']:
            molinit['units'] = units
        else:
            raise ValidationError('Invalid molecule geometry units: {}'.format(units))

        if input_units_to_au is not None:
            molinit['input_units_to_au'] = input_units_to_au

        if fix_com in [True, False]:
            molinit['fix_com'] = fix_com
        else:
            raise ValidationError('Invalid fix_com: {}'.format(fix_com))

        if fix_orientation in [True, False]:
            molinit['fix_orientation'] = fix_orientation
        else:
            raise ValidationError('Invalid fix_orientation: {}'.format(fix_com))

        if fix_symmetry is not None:
            molinit['fix_symmetry'] = fix_symmetry.lower()

        return molinit


def validate_and_fill_geometry(geom=None, tooclose=0.1):
    """Check `geom` for overlapping atoms. Return flattened"""

    if geom is None:
        raise ValidationError("""Geometry must be provided.""")

    npgeom = np.array(geom, dtype=np.float).reshape((-1, 3))
    dm = distance_matrix(npgeom, npgeom)

    iu = np.triu_indices(dm.shape[0])
    dm[iu] = 10.
    tooclosem = np.where(dm < tooclose)

    if tooclosem[0].shape[0]:
        raise ValidationError(
            """Following atoms are too close: {}""".format([(i, j, dm[i, j]) for i, j in zip(*tooclosem)]))

    return {'geom': npgeom.reshape((-1))}


def validate_and_fill_nuclei(nat,
                             elea=None,
                             elez=None,
                             elem=None,
                             mass=None,
                             real=None,
                             elbl=None,
                             # processing details
                             speclabel=True,
                             nonphysical=False,
                             mtol=1.e-3,
                             verbose=1):
    """Check the nuclear identity arrays for consistency and fill in knowable values."""

    if elea is None:
        elea = np.full((nat, ), None)
    else:
        # -1 equivalent to None
        elea = np.array([(None if at == -1 else at) for at in elea])

    if elez is None:
        elez = np.full((nat, ), None)
    else:
        elez = np.array(elez)

    if elem is None:
        elem = np.full((nat, ), None)
    else:
        elem = np.array(elem)

    if mass is None:
        mass = np.full((nat, ), None)
    else:
        mass = np.array(mass)

    if real is None:
        real = np.full((nat, ), None)
    else:
        real = np.array(real)

    if elbl is None:
        elbl = np.full((nat, ), None)
    else:
        elbl = np.array(elbl)

    if not ((nat, ) == elea.shape == elez.shape == elem.shape == mass.shape == real.shape == elbl.shape):
        raise ValidationError("""Dimension mismatch ({}) among A ({}), Z ({}), E ({}), mass ({}), real ({}), and elbl({})""".format(
            (nat,), elea.shape, elez.shape, elem.shape, mass.shape, real.shape, elbl.shape))

    A, Z, E, mass, real, label = zip(*[reconcile_nucleus(A=elea[at],
                                                         Z=elez[at],
                                                         E=elem[at],
                                                         mass=mass[at],
                                                         real=real[at],
                                                         label=elbl[at],
                                                         speclabel=speclabel,
                                                         nonphysical=nonphysical,
                                                         mtol=mtol,
                                                         verbose=verbose) for at in range(nat)])
    return {'elea': np.array(A, dtype=np.int),
            'elez': np.array(Z, dtype=np.int),
            'elem': np.array(E),
            'mass': np.array(mass, dtype=np.float),
            'real': np.array(real, dtype=np.bool),
            'elbl': np.array(label)}


def validate_and_fill_fragments(nat,
                                fragment_separators=None,
                                fragment_types=None,
                                fragment_charges=None,
                                fragment_multiplicities=None):
    """Check consistency of fragment specifiers wrt type and length. For
    charge & multiplicity, scientific defaults are not computed or applied;
    rather, missing slots are filled with `None` for later processing.

    """
    if fragment_separators is None:
        if fragment_types is None and fragment_charges is None and fragment_multiplicities is None:
            frs = [] #np.array([], dtype=np.int)  # if empty, needs to be both ndarray and int
            frt = ['Real']
            frc = [None]
            frm = [None]
        else:
            raise ValidationError(
                """Fragment quantities given without separation info: sep ({}), types ({}), chg ({}), and mult ({})""".
                format(fragment_separators, fragment_types, fragment_charges, fragment_multiplicities))
    else:
        trial_geom = np.zeros((nat, 3))
        try:
            split_geom = np.split(trial_geom, fragment_separators, axis=0)
        except TypeError:
            raise ValidationError("""fragment_separators ({}) unable to perform trial np.split on geometry.""".format(
                fragment_separators))
        if any(len(f) == 0 for f in split_geom):
            raise ValidationError(
                """fragment_separators ({}) yields zero-length fragment(s) after trial np.split on geometry.""".format(
                    split_geom))
        if sum(len(f) for f in split_geom) != nat:
            raise ValidationError(
                """fragment_separators ({}) yields overlapping fragment(s) after trial np.split on geometry, possibly unsorted.""".
                format(split_geom))
        frs = fragment_separators  #np.array(fragment_separators)
        nfr = len(split_geom)

        if fragment_types is None:
            frt = ['Real'] * nfr
        elif all(f in ['Real', 'Ghost', 'Absent'] for f in fragment_types):
            frt = fragment_types
        else:
            raise ValidationError("""fragment_types not among 'Real', 'Ghost', 'Absent': {}""".format(fragment_types))

        if fragment_charges is None:
            frc = [None] * nfr
        else:
            try:
                frc = [(f if f is None else float(f)) for f in fragment_charges]
            except TypeError:
                raise ValidationError("""fragment_charges not among None or float: {}""".format(fragment_charges))

        if fragment_multiplicities is None:
            frm = [None] * nfr
        elif all(f is None or (isinstance(f, (int, np.int64, long)) and f >= 1) for f in fragment_multiplicities):
            frm = fragment_multiplicities
        else:
            raise ValidationError(
                """fragment_multiplicities not among None or positive integer: {}""".format(fragment_multiplicities))

    if not (len(frt) == len(frc) == len(frm) == len(frs) + 1):
        raise ValidationError(
            """Dimension mismatch among fragment quantities: sep + 1 ({}), types ({}), chg ({}), and mult({})""".
            format(len(frs) + 1, len(frt), len(frc), len(frm)))

    return {'fragment_separators': list(frs),
            'fragment_charges': frc,
            'fragment_multiplicities': frm}


def from_string(molstr, dtype='psi4', return_processed=False):
    """

    Parameters
    ----------
    molstr : str
        Multiline string specification of molecule in a recognized format.
    dtype : {'xyz', 'xyz+', 'psi4'}, optional
        Molecule format name.
    return_processed : bool, optional
        Additionally return intermediate dictionary.

    Returns
    -------
    molrec : dict
        Molecule dictionary spec. See :py:func:`from_arrays`.
    molinit : dict, optional
        Intermediate "molrec"-like dictionary containing `molstr` info after
        processing by this function but before the validation and defaulting of
        `from_arrays` that returns the proper `molrec`.
        Only provided if `return_processed` is True.

    Raises
    ------
    ValidationError
        After processing of `molstr`, only an empty string should remain. Anything left is a syntax error.

    """
    ##"""Module with utility functions that act on molecule objects."""
    ##from psi4.driver.inputparser import process_pubchem_command, pubchemre
    #molecule = init_psi4_molecule_from_any_string(geom, name=name)

    ##    #if not has_efp:
    #    # Figure out how and if we will parse the Molecule adata
    #    mname = kwargs.pop("name", "default")
    #    dtype = kwargs.pop("dtype", "psi4").lower()
    #    if mol_str is not None:
    # << 2 >> dict[m] -- seed Psi4 minimal and default fields

    print('<<< QMMOL', molstr, '>>>')

    # << 1 >> str -- discard comments
    molstr = filter_comments(molstr)

    if dtype == 'xyz':
        """Strict XYZ file format

        Specifiable: geom, elem/elez (element identity)
        Inaccessible: mass, real (vs. ghost), elbl (user label), name, units (assumed [A]), input_units_to_au, fix_com/orientation/symmetry, fragmentation, molecular_charge, molecular_multiplicity

        <number of atoms>
        comment line
        <element_symbol or atomic_number> <x> <y> <z>
        ...
        <element_symbol or atomic_number> <x> <y> <z>

        """

    elif dtype == 'xyz+':
        """Enhanced XYZ format

        <number of atoms> [<[bohr|au|ang]>]
        [<molecular_charge> <molecular_multiplicity>] comment line
        <psi4_nucleus_spec> <x> <y> <z>
        ...
        <psi4_nucleus_spec> <x> <y> <z>

        Specifiable: geom, elem/elez (element identity), mass, real (vs. ghost), elbl (user label), units (defaults [A]), molecular_charge, molecular_multiplicity
        Inaccessible: name, input_units_to_au, fix_com/orientation/symmetry, fragmentation

        """

    elif dtype == 'psi4':

        molinit = {}

        # << 2 >>  str-->dict[q] -- process units, com, orient, symm
        molstr, processed = _filter_universals(molstr)
        molinit.update(processed)

        #        # << 3 >>  str-->dict[e] -- process efp
        #        mol_str, efp_init = filter_libefp(mol_str, confer=mol_init)
        #        if efp_init:
        #            print('<<< core.EFP INTO', efp_init, '>>>')
        #            # GOOD! core.efp_init()
        #            # GOOD! efp = core.get_active_efp()
        #            # GOOD! efp.construct_from_pydict(efp_init)
        #            # << 4 >>  dict[e]-->dict[m] --- tie qm & efp axes
        #            mol_init['fix_com'] = True
        #            mol_init['fix_orientation'] = True
        #            mol_init['fix_symmetry'] = 'c1'

        # << 3 >> str-->dict[q] -- process atoms, chg, mult, frags
        molstr, processed = _filter_mints(molstr)
        molinit.update(processed)

        if molstr:
            print('\n<<< IFNAL:\n', molstr, '\n>>>\n')
            raise ValidationError("""Unprocessable Molecule remanents:\n\t{}""".format(molstr))

        import pprint
        pprint.pprint(molinit)

        print('\nINTO from_arrays <<<', molinit, '>>>\n')
        molrec = from_arrays(**molinit, speclabel=True)
        print('\nMOLREC <<<', molrec, '>>>\n')

        #        return {'molecule': mol_init, 'libefp': efp_init}
        #        #has_efp, geom = filter_libefp(geom)

        if return_processed:
            return molrec, molinit
        else:
            return molrec

    else:
        raise KeyError("Molecule: dtype of %s not recognized.")


#def filter_pubchem(mol_str):
#    pass
#    pubchemerror = re.compile(r'^\s*PubchemError\s*$', re.IGNORECASE)
#    pubcheminput = re.compile(r'^\s*PubchemInput\s*$', re.IGNORECASE)


def _filter_universals(string):
    """Process multiline `string` for fix_ and unit markers,
    returning a string of unprocessed `string` and a dictionary of
    processed fields.

    fix_com
    fix_orientation
    fix_symmetry
#    input_units_to_au (not settable)
    units

    """
    com = re.compile(r'\A(no_com|nocom)\Z', re.IGNORECASE)
    orient = re.compile(r'\A(no_reorient|noreorient)\Z', re.IGNORECASE)
    bohrang = re.compile(r'\Aunits?[\s=]+((?P<ubohr>(bohr|au|a.u.))|(?P<uang>(ang|angstrom)))\Z', re.IGNORECASE)
    symmetry = re.compile(r'\Asymmetry[\s=]+(?P<pg>\w+)\Z', re.IGNORECASE)

    def process_com(matchobj):
        processed['fix_com'] = True
        return ''

    def process_orient(matchobj):
        processed['fix_orientation'] = True
        return ''

    def process_bohrang(matchobj):
        if matchobj.group('uang'):
            processed['units'] = 'Angstrom'
        elif matchobj.group('ubohr'):
            processed['units'] = 'Bohr'
        return ''

    def process_symmetry(matchobj):
        processed['fix_symmetry'] = matchobj.group('pg').lower()
        return ''

    reconstitute = []
    processed = {}

    for line in string.split('\n'):
        line = re.sub(com, process_com, line.strip())
        line = re.sub(orient, process_orient, line)
        line = re.sub(bohrang, process_bohrang, line)
        line = re.sub(symmetry, process_symmetry, line)
        if line:
            reconstitute.append(line)

    return '\n'.join(reconstitute), processed


def _filter_mints(string):
    """Handle extracting fragment, atom, and chg/mult lines from `string`.

    Returns
    -------
    str, dict
        Returns first a subset (plus some fragment separation guidance) of
            `string` containing the unmatched contents. These are generally input
            violations unless handled by a subsequent processing function.
        Returns second a dictionary with processed extractions. Contains (some
            optional) the following keys.

            molecular_charge : float, optional
            molecular_multiplicity : int, optional
            geom
            elbl
            fragment_separators
            fragment_charges
            fragment_multiplicities

    """
    fragment_marker = re.compile(r'^\s*--\s*$', re.MULTILINE)
    CHG = r'(?P<chg>' + NUMBER + r')'
    MULT = r'(?P<mult>\d+)'
    cgmp = re.compile(r'\A' + CHG + SEP + MULT + r'\Z', re.VERBOSE)
    atom_cartesian = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + r'(?P<x>' + NUMBER + r')' + 
                                                                           SEP + r'(?P<y>' + NUMBER + r')' +
                                                                           SEP + r'(?P<z>' + NUMBER + r')' + r'\Z', re.IGNORECASE | re.VERBOSE)
    atom_zmat1 = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + r'\Z', re.IGNORECASE | re.VERBOSE)
    atom_zmat2 = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + r'(?P<Ridx>\d+)' + SEP + r'(?P<Rval>' + NUMBER + r')' + r'\Z', re.IGNORECASE | re.VERBOSE)
    atom_zmat3 = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + r'(?P<Ridx>\d+)' + SEP + r'(?P<Rval>' + NUMBER + r')' + 
                                                                       SEP + r'(?P<Aidx>\d+)' + SEP + r'(?P<Aval>' + NUMBER + r')' + r'\Z', re.IGNORECASE | re.VERBOSE)
    atom_zmat4 = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + r'(?P<Ridx>\d+)' + SEP + r'(?P<Rval>' + NUMBER + r')' + 
                                                                       SEP + r'(?P<Aidx>\d+)' + SEP + r'(?P<Aval>' + NUMBER + r')' +
                                                                       SEP + r'(?P<Didx>\d+)' + SEP + r'(?P<Dval>' + NUMBER + r')' + r'\Z', re.IGNORECASE | re.VERBOSE)
    #variable = re.compile(r'^\s*(\w+)\s*=\s*(-?\d+\.\d+|-?\d+\.|-?\.\d+|-?\d+|tda)\s*$', re.IGNORECASE)


    def process_system_cgmp(matchobj):
        """Handles optional special first fragment with sole contents overall chg/mult."""

        processed['molecular_charge'] = float(matchobj.group('chg'))
        processed['molecular_multiplicity'] = int(matchobj.group('mult'))
        return ''

    def filter_fragment(fstring):
        """Handles extraction from everything within a fragment marker "--" of a
        single chg/mult (or None/None) and multiple atom lines.

        """

        def process_fragment_cgmp(matchobj):
            print('frcgmp hit', matchobj.groups())
            processed['fragment_charges'].append(float(matchobj.group('chg')))
            processed['fragment_multiplicities'].append(int(matchobj.group('mult')))
            return ''

        def process_atom_cartesian(matchobj):
            #atom_init['qm_type'] = 'qmcart'
            #atom_init['z'] = float(matchobj.group(10)) * input_units_to_au

            processed['elbl'].append(matchobj.group('nucleus'))
            processed['geom'].append(float(matchobj.group('x')))
            processed['geom'].append(float(matchobj.group('y')))
            processed['geom'].append(float(matchobj.group('z')))
            return ''

        freconstitute = []
        start_atom = len(processed["elbl"])
        if start_atom > 0:
            processed['fragment_separators'].append(start_atom)

        fcgmp_found = False
        for iln, line in enumerate(fstring.split('\n')):
            if not fcgmp_found:
                line, fcgmp_found = re.subn(cgmp, process_fragment_cgmp, line.strip())
            line = re.sub(atom_cartesian, process_atom_cartesian, line)
            #line = re.sub(atom_zmat1, process_atom_zmat1, line)
            if line:
                freconstitute.append(line)

        if not fcgmp_found:
            processed['fragment_charges'].append(None)
            processed['fragment_multiplicities'].append(None)

        return '\n'.join(freconstitute), processed

    reconstitute = []
    processed = {}
    processed['geom'] = []
    processed['elbl'] = []
    processed['fragment_separators'] = []
    processed['fragment_charges'] = []
    processed['fragment_multiplicities'] = []

    # handle `--`-demarcated blocks
    for ifr, frag in enumerate(re.split(fragment_marker, string)):
        frag = frag.strip()
        if ifr == 0 and cgmp.match(frag):
            frag, ntotch = re.subn(cgmp, process_system_cgmp, frag)
        else:
            frag, processed = filter_fragment(frag)
        if frag:
            reconstitute.append(frag)

    return '\n-- (guidance)\n'.join(reconstitute), processed
