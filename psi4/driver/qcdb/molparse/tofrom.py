import re
import pprint

import numpy as np

from ..util import distance_matrix, update_with_error
from ..exceptions import *
from ..physconst import psi_bohr2angstroms
from .chgmult import validate_and_fill_chgmult
from .nucleus import reconcile_nucleus
from .regex import *
from . import pubchem

try:
    long(1)
except NameError:
    long = int


def from_input_arrays(enable_qm=True,
                      enable_efp=True,
                      missing_enabled_return_qm='error',
                      missing_enabled_return_efp='error',
                      # qm
                      geom=None,
                      elea=None,
                      elez=None,
                      elem=None,
                      mass=None,
                      real=None,
                      elbl=None,
                      name=None,
                      units='Angstrom',
                      input_units_to_au=None,
                      fix_com=None,
                      fix_orientation=None,
                      fix_symmetry=None,
                      fragment_separators=None,
                      fragment_charges=None,
                      fragment_multiplicities=None,
                      molecular_charge=None,
                      molecular_multiplicity=None,
                      # efp
                      fragment_files=None,
                      hint_types=None,
                      geom_hints=None,
                      # processing details
                      speclabel=True,
                      tooclose=0.1,
                      zero_ghost_fragments=False,
                      nonphysical=False,
                      mtol=1.e-3,
                      verbose=1):

    molinit = {}
    if enable_qm:
        molinit['qm'] = {}
    if enable_efp:
        molinit['efp'] = {}

    if enable_efp:
        processed = from_arrays(domain='efp',
                                missing_enabled_return=missing_enabled_return_efp,
                                units=units,
                                input_units_to_au=input_units_to_au,
                                fix_com=fix_com,
                                fix_orientation=fix_orientation,
                                fix_symmetry=fix_symmetry,
                                fragment_files=fragment_files,
                                hint_types=hint_types,
                                geom_hints=geom_hints,
                                # which other processing details needed?
                                verbose=verbose)
        update_with_error(molinit, {'efp': processed})
        if molinit['efp'] == {}:
            del molinit['efp']

    efp_present = enable_efp and 'efp' in molinit and bool(len(molinit['efp']['geom_hints']))
    if efp_present:
        fix_com = True
        fix_orientation = True
        fix_symmetry = 'c1'

    if enable_qm:
        processed = from_arrays(domain='qm',
                                missing_enabled_return=missing_enabled_return_qm,
                                geom=geom,
                                elea=elea,
                                elez=elez,
                                elem=elem,
                                mass=mass,
                                real=real,
                                elbl=elbl,
                                name=name,
                                units=units,
                                input_units_to_au=input_units_to_au,
                                fix_com=fix_com,
                                fix_orientation=fix_orientation,
                                fix_symmetry=fix_symmetry,
                                fragment_separators=fragment_separators,
                                fragment_charges=fragment_charges,
                                fragment_multiplicities=fragment_multiplicities,
                                molecular_charge=molecular_charge,
                                molecular_multiplicity=molecular_multiplicity,
                                # processing details
                                speclabel=speclabel,
                                tooclose=tooclose,
                                zero_ghost_fragments=zero_ghost_fragments,
                                nonphysical=nonphysical,
                                mtol=mtol,
                                verbose=1)
        update_with_error(molinit, {'qm': processed})
        if molinit['qm'] == {}:
            del molinit['qm']

    return molinit


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
                fix_com=None,
                fix_orientation=None,
                fix_symmetry=None,
                fragment_separators=None,
                fragment_charges=None,
                fragment_multiplicities=None,
                molecular_charge=None,
                molecular_multiplicity=None,
                fragment_files=None,
                hint_types=None,
                geom_hints=None,

                domain='qm',
                missing_enabled_return='error',

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
    missing_enabled_return : {'minimal', 'none', 'error'}
        What to do when an enabled domain is of zero-length? Respectively, return
        a fully valid but empty molrec, return empty dictionary, or throw error.


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

        (5) apart from some mild optional fields, _all_ fields will
        be present (correlary of "fully validated and defaulted") - no
        need to check for every key. in some cases like efp, keys will
        appear in blocks, so pre-handshake there will be a few hint keys
        and post-handshake they will be joined by full qm-like molrec.

        (6) molrec should be idempotent through this function (equiv to
        schema validator) but are not idempostent throughout its life. if
        fields permit, frame may be changed. Future? if fields permit,
        mol may be symmetrized. Coordinates and angles may change units
        or range if program returns them in only one form.

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

    EFP extension (this + units is minimal)

    fragment_files : list of str
        (nfr, ) lowercased names of efp meat fragment files.
    hint_types : {'xyzabc', 'points'}
        (nfr, ) type of fragment orientation hint.
    geom_hints : list of lists of float
        (nfr, ) inner lists have length 6 (xyzabc; to orient the center) or
        9 (points; to orient the first three atoms) of the EFP fragment.

    """
    # <<  domain sorting  >>
    available_domains = ['qm', 'efp']
    if domain not in available_domains:
        raise ValidationError('Topology domain {} not available for processing. Choose among {}'.format(domain, available_domains))

    if domain == 'qm' and geom is None or geom == []:
        if missing_enabled_return == 'none':
            return {}
        elif missing_enabled_return == 'minimal':
            geom = []
        else:
            raise ValidationError("""For domain 'qm', `geom` must be provided.""")
    if domain == 'efp' and geom_hints is None or geom_hints == []:
        if missing_enabled_return == 'none':
            return {}
        elif missing_enabled_return == 'minimal':
            geom_hints = []
            fragment_files = []
            hint_types = []
        else:
            raise ValidationError("""For domain 'efp', `geom_hints` must be provided.""")

    molinit = {}
    extern = False

    processed = validate_and_fill_units(name=name,
                                        units=units,
                                        input_units_to_au=input_units_to_au,
                                        always_return_iutau=False)
    update_with_error(molinit, processed)

    if domain == 'efp':
        processed = validate_and_fill_efp(fragment_files=fragment_files,
                                          hint_types=hint_types,
                                          geom_hints=geom_hints)
        update_with_error(molinit, processed)
        extern = bool(len(molinit['geom_hints']))

    if domain == 'qm' or (domain == 'efp' and geom is not None):
        processed = validate_and_fill_geometry(geom=geom,
                                               tooclose=tooclose)
        update_with_error(molinit, processed)
        nat = molinit['geom'].shape[0] // 3

        processed = validate_and_fill_nuclei(nat,
                                             elea=elea,
                                             elez=elez,
                                             elem=elem,
                                             mass=mass,
                                             real=real,
                                             elbl=elbl,
                                             speclabel=speclabel,
                                             nonphysical=nonphysical,
                                             mtol=mtol,
                                             verbose=verbose)
        update_with_error(molinit, processed)

        processed = validate_and_fill_fragments(nat,
                                                fragment_separators=fragment_separators,
                                                fragment_charges=fragment_charges,
                                                fragment_multiplicities=fragment_multiplicities)
        update_with_error(molinit, processed)

        Z_available = molinit['elez'] * molinit['real'] * 1.
        processed = validate_and_fill_chgmult(zeff=Z_available,
                                              fragment_separators=molinit['fragment_separators'],
                                              molecular_charge=molecular_charge,
                                              fragment_charges=molinit['fragment_charges'],
                                              molecular_multiplicity=molecular_multiplicity,
                                              fragment_multiplicities=molinit['fragment_multiplicities'],
                                              zero_ghost_fragments=zero_ghost_fragments,
                                              verbose=verbose)
        del molinit['fragment_charges']  # sometimes safe update is too picky about overwriting v_a_f_fragments values
        del molinit['fragment_multiplicities']
        update_with_error(molinit, processed)

    extern = (domain == 'efp')

    processed = validate_and_fill_frame(extern=extern,
                                        fix_com=fix_com,
                                        fix_orientation=fix_orientation,
                                        fix_symmetry=fix_symmetry)
    update_with_error(molinit, processed)

    print('RETURN FROM from_arrays', domain.upper())
    pprint.pprint(molinit)

    return molinit



def validate_and_fill_units(name=None,
                            units='Angstrom',
                            input_units_to_au=None,
                            always_return_iutau=False):
    molinit = {}

    if name is not None:
        molinit['name'] = name

    if units in ['Angstrom', 'Bohr']:
        molinit['units'] = units
    else:
        raise ValidationError('Invalid molecule geometry units: {}'.format(units))

    if molinit['units'] == 'Bohr':
        iutau = 1.
    elif molinit['units'] == 'Angstrom':
        iutau = 1. / psi_bohr2angstroms

    if input_units_to_au is not None:
        if abs(input_units_to_au - iutau) < 0.05:
            iutau = input_units_to_au
        else:
            raise ValidationError("""No big perturbations to physical constants! {} !~= {}""".format(iutau, input_units_to_au))

    if always_return_iutau or input_units_to_au is not None:
        molinit['input_units_to_au'] = iutau

    return molinit


def validate_and_fill_frame(extern,
                            fix_com=None,
                            fix_orientation=None,
                            fix_symmetry=None):

    if fix_com is True:
        com = True
    elif fix_com is False:
        if extern:
            raise ValidationError('Invalid fix_com ({}) with extern ({})'.format(fix_com, extern))
        else:
            com = False
    elif fix_com is None:
        com = extern
    else:
        raise ValidationError('Invalid fix_com: {}'.format(fix_com))

    if fix_orientation is True:
        orient = True
    elif fix_orientation is False:
        if extern:
            raise ValidationError('Invalid fix_orientation ({}) with extern ({})'.format(fix_orientation, extern))
        else:
            orient = False
    elif fix_orientation is None:
        orient = extern
    else:
        raise ValidationError('Invalid fix_orientation: {}'.format(fix_orientation))

    symm = None
    if extern:
        if fix_symmetry is None:
            symm = 'c1'
        elif fix_symmetry.lower() == 'c1':
            symm = 'c1'
        else:
            raise ValidationError('Invalid (non-C1) fix_symmetry ({}) with extern ({})'.format(fix_symmetry, extern))
    else:
        if fix_symmetry is not None:
            symm = fix_symmetry.lower()

    molinit = {}
    molinit['fix_com'] = com
    molinit['fix_orientation'] = orient
    if symm:
        molinit['fix_symmetry'] = symm

    return molinit


def preserve_validate_and_fill_frame(extern,
                            fix_com=None,
                            fix_orientation=None,
                            fix_symmetry=None):

    if fix_com is True:
        com = True
    elif fix_com is False:
        if extern:
            raise ValidationError('Invalid fix_com ({}) with extern ({})'.format(fix_com, extern))
        else:
            com = False
    elif fix_com is None:
        com = extern
    else:
        raise ValidationError('Invalid fix_com: {}'.format(fix_com))

    if fix_orientation is True:
        orient = True
    elif fix_orientation is False:
        if extern:
            raise ValidationError('Invalid fix_orientation ({}) with extern ({})'.format(fix_orientation, extern))
        else:
            orient = False
    elif fix_orientation is None:
        orient = extern
    else:
        raise ValidationError('Invalid fix_orientation: {}'.format(fix_orientation))

    symm = None
    if extern:
        if fix_symmetry is None:
            symm = 'c1'
        elif fix_symmetry.lower() == 'c1':
            symm = 'c1'
        else:
            raise ValidationError('Invalid (non-C1) fix_symmetry ({}) with extern ({})'.format(fix_symmetry, extern))
    else:
        if fix_symmetry is not None:
            symm = fix_symmetry.lower()

    molinit = {}
    molinit['fix_com'] = com
    molinit['fix_orientation'] = orient
    if symm:
        molinit['fix_symmetry'] = symm

    return molinit


def validate_and_fill_efp(fragment_files=None,
                          hint_types=None,
                          geom_hints=None):
    """

    Parameters
    ----------
#    standardize_efp_angles_units : bool, optional
#        Move abc of xyzabc hints into (-pi, pi] range returned by libefp.
#        Not needed for input as libefp takes any range as input but helpful for
#        matching output.

    """
    if (fragment_files is None or hint_types is None or geom_hints is None or
        fragment_files == [None] or hint_types == [None] or geom_hints == [None] or
        not (len(fragment_files) == len(hint_types) == len(geom_hints))):

        raise ValidationError(
                """Missing or inconsistent length among efp quantities: fragment_files ({}), hint_types ({}), and geom_hints ({})""".
                format(fragment_files, hint_types, geom_hints))

    # NOTE: imposing case on file
    try:
        files = [f.lower() for f in fragment_files]
    except AttributeError:
        raise ValidationError("""fragment_files not strings: {}""".format(fragment_files))

    if all(f in ['xyzabc', 'points', 'rotmat'] for f in hint_types):
        types = hint_types
    else:
        raise ValidationError("""hint_types not among 'xyzabc', 'points', 'rotmat': {}""".format(hint_types))

    hints = []
    hlen = {'xyzabc': 6,
            'points': 9,
            'rotmat': 12}
    for ifr, fr in enumerate(geom_hints):
        try:
            hint = [float(f) for f in fr]
        except ValueError:
            raise ValidationError("""Un float-able elements in geom_hints[{}]: {}""".format(ifr, fr))

        htype = hint_types[ifr]
        if len(hint) == hlen[htype]:
            hints.append(hint)
        else:
            raise ValidationError("""EFP hint type {} not {} elements: {}""".format(
                htype, hlen[htype], hint))

    return {'fragment_files': files,
            'hint_types': types,
            'geom_hints': hints}


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

    if nat:
        A, Z, E, mass, real, label = zip(*[reconcile_nucleus(A=elea[at],
    #theans = zip(*[reconcile_nucleus(A=elea[at],
                                                         Z=elez[at],
                                                         E=elem[at],
                                                         mass=mass[at],
                                                         real=real[at],
                                                         label=elbl[at],
                                                         speclabel=speclabel,
                                                         nonphysical=nonphysical,
                                                         mtol=mtol,
                                                         verbose=verbose) for at in range(nat)])
    else:
        A = Z = E = mass = real = label = []
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
            if nat != 0:
                raise ValidationError(
                    """fragment_separators ({}) yields zero-length fragment(s) after trial np.split on geometry.""".format(
                        split_geom))
        if sum(len(f) for f in split_geom) != nat:
            raise ValidationError(
                """fragment_separators ({}) yields overlapping fragment(s) after trial np.split on geometry, possibly unsorted.""".
                format(split_geom))
        frs = fragment_separators
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

    #if nat == 0:
    #    return {'fragment_separators': [],
    #            'fragment_charges': [],
    #            'fragment_multiplicities': []}

    if not (len(frt) == len(frc) == len(frm) == len(frs) + 1):
        raise ValidationError(
            """Dimension mismatch among fragment quantities: sep + 1 ({}), types ({}), chg ({}), and mult({})""".
            format(len(frs) + 1, len(frt), len(frc), len(frm)))

    return {'fragment_separators': list(frs),
            'fragment_charges': frc,
            'fragment_multiplicities': frm}


def from_string(molstr, dtype='psi4',
                return_processed=False,
                enable_qm=True, 
                enable_efp=True,
                missing_enabled_return_qm='none',
                missing_enabled_return_efp='none'):
    """

    Parameters
    ----------
    molstr : str
        Multiline string specification of molecule in a recognized format.
    dtype : {'xyz', 'xyz+', 'psi4'}, optional
        Molecule format name.
    return_processed : bool, optional
        Additionally return intermediate dictionary.
    enable_qm : bool, optional
    enable_efp: bool, optional

    Returns
    -------
    molrec : dict
        Molecule dictionary spec. See :py:func:`from_arrays`.
    molinit : dict, optional
        Intermediate "molrec"-like dictionary containing `molstr` info after
        parsing by this function but before the validation and defaulting of
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

    #    # Figure out how and if we will parse the Molecule adata
    #    mname = kwargs.pop("name", "default")
    #    dtype = kwargs.pop("dtype", "psi4").lower()

    molinit = {}

    print('<<< QMMOL', molstr, '>>>')

    # << 1 >>  str-->str -- discard comments
    molstr = filter_comments(molstr)

    if dtype == 'xyz':
        """Strict XYZ file format

        String Layout
        -------------
        <number of atoms>
        comment line
        <element_symbol or atomic_number> <x> <y> <z>
        ...
        <element_symbol or atomic_number> <x> <y> <z>

        QM Domain
        ---------
        Specifiable: geom, elem/elez (element identity)
        Inaccessible: mass, real (vs. ghost), elbl (user label), name, units (assumed [A]),
                      input_units_to_au, fix_com/orientation/symmetry, fragmentation,
                      molecular_charge, molecular_multiplicity

        """
        # << 2 >>  str-->dict -- process atoms, units
        molstr, processed = _filter_xyz(molstr, strict=True)
        molinit.update(processed)

    elif dtype == 'xyz+':
        """Enhanced XYZ format

        String Layout
        -------------
        <number of atoms> [<bohr|au|ang>]
        [<molecular_charge> <molecular_multiplicity>] comment line
        <psi4_nucleus_spec> <x> <y> <z>
        ...
        <psi4_nucleus_spec> <x> <y> <z>

        QM Domain
        ---------
        Specifiable: geom, elem/elez (element identity), mass, real (vs. ghost), elbl (user label),
                     units (defaults [A]), molecular_charge, molecular_multiplicity
        Inaccessible: name, input_units_to_au, fix_com/orientation/symmetry, fragmentation

        Note that <number of atoms> is pattern-matched but ignored.

        """
        # << 2 >>  str-->dict -- process atoms, units, chg, mult
        molstr, processed = _filter_xyz(molstr, strict=False)
        molinit.update(processed)

    elif dtype == 'psi4':
        """Psi4 molecule {...} format

        QM Domain
        ---------
        Specifiable: geom, elem/elez (element identity), mass, real (vs. ghost), elbl (user label),
                     units (defaults [A]), fix_com/orientation/symmetry, fragment_separators,
                     fragment_charges, fragment_multiplicities, molecular_charge, molecular_multiplicity
        Inaccessible: name, input_units_to_au

        EFP Domain
        ----------
        Specifiable: units, fix_com/orientation/symmetry, fragment_files, hint_types, geom_hints
        Inaccessible: anything atomic or fragment details -- geom, elem/elez (element identity),
                      mass, real (vs. ghost), elbl (user label), fragment_separators, fragment_charges,
                      fragment_multiplicities, molecular_charge, molecular_multiplicity

        """
        # Notes
        # * *_filter functions must fill non-overlapping fields
        # * not recc but can add to downstream by appending to str

        # << 2.1 >>  str-->str -- process pubchem into str for downfunction
        molstr, processed = _filter_pubchem(molstr)
        molinit.update(processed)

        # << 2.2 >>  str-->dict -- process units, com, orient, symm
        molstr, processed = _filter_universals(molstr)
        molinit.update(processed)

        # << 2.3 >>  str-->dict -- process efp frags
        molstr, processed = _filter_libefp(molstr)
        molinit.update(processed)

        # << 2.4 >>  str-->dict -- process atoms, chg, mult, frags
        molstr, processed = _filter_mints(molstr)
        molinit.update(processed)

    else:
        raise KeyError("Molecule: dtype of %s not recognized.")

    if molstr:
        raise ValidationError("""Unprocessable Molecule remanents:\n\t{}""".format(molstr))

    print('\nINTO from_input_arrays <<<')
    pprint.pprint(molinit)
    print('>>>\n')

    # << 3 >>  dict-->molspec
    molrec = from_input_arrays(**molinit,
                               speclabel=True,
                               enable_qm=enable_qm, 
                               enable_efp=enable_efp,
                               missing_enabled_return_qm=missing_enabled_return_qm,
                               missing_enabled_return_efp=missing_enabled_return_efp)

    print('\nMOLREC <<<', molrec, '>>>\n')

    if return_processed:
        return molrec, molinit
    else:
        return molrec


#    pubchemerror = re.compile(r'^\s*PubchemError\s*$', re.IGNORECASE)
#    pubcheminput = re.compile(r'^\s*PubchemInput\s*$', re.IGNORECASE)
#        # N.B. Anything starting with PubchemError will be handled correctly by the molecule parser
#        # in libmints, which will just print the rest of the string and exit gracefully.

def _filter_pubchem(string):
    pubchemre = re.compile(r'\Apubchem' + r'\s*:\s*' + r'(?P<pubsearch>(\S+))\Z', re.IGNORECASE)

    def process_pubchem(matchobj):
        pubsearch = matchobj.group('pubsearch')

        if pubsearch.isdigit():
            # just a number - must be a CID
            pcobj = pubchem.PubChemObj(int(pubsearch), '', '')
            try:
                xyz = pcobj.getMoleculeString()
            except Exception as e:
                raise ValidationError(e.message)
        else:
            # search pubchem for the provided string
            try:
                results = pubchem.getPubChemResults(pubsearch)
            except Exception as e:
                raise ValidationError(e.message)

            if not results:
                # Nothing!
                raise ValidationError("""PubchemError: No results were found when searching PubChem for {}.""".format(pubsearch))
            elif len(results) == 1:
                # There's only 1 result - use it
                xyz = results[0].getMoleculeString()
            else:
                # There are multiple results
                for result in results:
                    if result.name().lower() == pubsearch.lower():
                        # We've found an exact match!
                        xyz = result.getMoleculeString()
                else:
                    # There are multiple results and none exact. Print and exit
                    msg = "\tPubchemError\n"
                    msg += "\tMultiple pubchem results were found. Replace\n\n\t\tpubchem:%s\n\n" % (pubsearch)
                    msg += "\twith the Chemical ID number or exact name from one of the following and re-run.\n\n"
                    msg += "\t Chemical ID     IUPAC Name\n\n"
                    for result in results:
                        msg += "%s" % (result)
                    raise ValidationErrror(msg)

        # remove PubchemInput first line and assert [A]
        xyz = xyz.replace('PubchemInput', 'units ang')
        return xyz

    reconstitute = []
    processed = {}

    for line in string.split('\n'):
        line = re.sub(pubchemre, process_pubchem, line.strip())
        if line:
            reconstitute.append(line)

    return '\n'.join(reconstitute), processed


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
    com_found = False
    orient_found = False
    bohrang_found = False
    symmetry_found = False

    for line in string.split('\n'):
        line = line.strip()
        if not com_found:
            line, com_found = re.subn(com, process_com, line)
        if not orient_found:
            line, orient_found = re.subn(orient, process_orient, line)
        if not bohrang_found:
            line, bohrang_found = re.subn(bohrang, process_bohrang, line)
        if not symmetry_found:
            line, symmetry_found = re.subn(symmetry, process_symmetry, line)
        if line:
            reconstitute.append(line)

    return '\n'.join(reconstitute), processed


def _filter_libefp(string):

    fragment_marker = re.compile(r'^\s*--\s*$', re.MULTILINE)
    efpxyzabc = re.compile(
        r'\A' + r'efp' + SEP + r'(?P<efpfile>(\w+))' + SEP +
        r'(?P<x>' + NUMBER + r')' + SEP + r'(?P<y>' + NUMBER + r')' + SEP + r'(?P<z>' + NUMBER + r')' + SEP +
        r'(?P<a>' + NUMBER + r')' + SEP + r'(?P<b>' + NUMBER + r')' + SEP + r'(?P<c>' + NUMBER + r')' + ENDL + r'\Z',
        re.IGNORECASE | re.VERBOSE)
    efppoints = re.compile(
        r'\A' + r'efp' + SEP + r'(?P<efpfile>(\w+))' + ENDL +
        r'[\s,]*' + r'(?P<x1>' + NUMBER + r')' + SEP + r'(?P<y1>' + NUMBER + r')' + SEP + r'(?P<z1>' + NUMBER + r')' + ENDL +
        r'[\s,]*' + r'(?P<x2>' + NUMBER + r')' + SEP + r'(?P<y2>' + NUMBER + r')' + SEP + r'(?P<z2>' + NUMBER + r')' + ENDL +
        r'[\s,]*' + r'(?P<x3>' + NUMBER + r')' + SEP + r'(?P<y3>' + NUMBER + r')' + SEP + r'(?P<z3>' + NUMBER + r')' + ENDL + r'\Z',
        re.IGNORECASE | re.MULTILINE | re.VERBOSE)

    def process_efpxyzabc(matchobj):
        processed['fragment_files'].append(matchobj.group('efpfile'))
        processed['hint_types'].append('xyzabc')
        processed['geom_hints'].append([
            float(matchobj.group('x')), float(matchobj.group('y')), float(matchobj.group('z')),
            float(matchobj.group('a')), float(matchobj.group('b')), float(matchobj.group('c'))])
        return ''

    def process_efppoints(matchobj):
        processed['fragment_files'].append(matchobj.group('efpfile'))
        processed['hint_types'].append('points')
        processed['geom_hints'].append([
            float(matchobj.group('x1')), float(matchobj.group('y1')), float(matchobj.group('z1')),
            float(matchobj.group('x2')), float(matchobj.group('y2')), float(matchobj.group('z2')),
            float(matchobj.group('x3')), float(matchobj.group('y3')), float(matchobj.group('z3'))])
        return ''

    reconstitute = []
    processed = {}
    processed['fragment_files'] = []
    processed['hint_types'] = []
    processed['geom_hints'] = []

    # handle `--`-demarcated blocks
    for frag in re.split(fragment_marker, string):
        frag = re.sub(efpxyzabc, process_efpxyzabc, frag.strip())
        frag = re.sub(efppoints, process_efppoints, frag)
        if frag:
            reconstitute.append(frag)

    return '\n--\n'.join(reconstitute), processed



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
    cgmp = re.compile(r'\A' + CHGMULT + r'\Z', re.VERBOSE)
    atom_cartesian = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + CARTXYZ + r'\Z', re.IGNORECASE | re.VERBOSE)
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
            processed['fragment_charges'].append(float(matchobj.group('chg')))
            processed['fragment_multiplicities'].append(int(matchobj.group('mult')))
            return ''

        def process_atom_cartesian(matchobj):
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

    return '\n--\n'.join(reconstitute), processed


def _filter_xyz(string, strict):
    """Handle extracting atom, units, and chg/mult lines from `string`.

    Parameters
    ----------
    strict : bool
        Whether to enforce a strict XYZ file format or to allow units, chg/mult,
        and add'l atom info.

    Returns
    -------
    str, dict
        Returns first a subset of `string` containing the unmatched contents.
        These are generally input violations unless handled by a subsequent
        processing function.
        Returns second a dictionary with processed extractions. Contains (some
            optional) the following keys.

            molecular_charge : float, optional (`strict=False` only)
            molecular_multiplicity : int, optional (`strict=False` only)
            geom
            elbl
            units : {'Angstrom', 'Bohr'} (`Bohr` `strict=False` only)

    """
    xyz1strict = re.compile(r'\A' + r'(?P<nat>\d+)' + r'\Z')
    SIMPLENUCLEUS = r"""((?P<E>[A-Z]{1,3})|(?P<Z>\d{1,3}))"""
    atom_cartesian_strict = re.compile(r'\A' + r'(?P<nucleus>' + SIMPLENUCLEUS + r')' + SEP + CARTXYZ + r'\Z', re.IGNORECASE | re.VERBOSE)

    xyz1 = re.compile(r'\A' + r'(?P<nat>\d+)' +
                      r'[\s,]*' + r'((?P<ubohr>(bohr|au))|(?P<uang>ang))?' + r'\Z', re.IGNORECASE)
    xyz2 = re.compile(r'\A' + CHGMULT, re.VERBOSE)
    atom_cartesian = re.compile(r'\A' + r'(?P<nucleus>' + NUCLEUS + r')' + SEP + CARTXYZ + r'\Z', re.IGNORECASE | re.VERBOSE)

    def process_bohrang(matchobj):
        nat = matchobj.group('nat')
        if matchobj.group('uang'):
            processed['units'] = 'Angstrom'
        elif matchobj.group('ubohr'):
            processed['units'] = 'Bohr'
        return ''

    def process_system_cgmp(matchobj):
        processed['molecular_charge'] = float(matchobj.group('chg'))
        processed['molecular_multiplicity'] = int(matchobj.group('mult'))
        return ''

    def process_atom_cartesian(matchobj):
        processed['elbl'].append(matchobj.group('nucleus'))
        processed['geom'].append(float(matchobj.group('x')))
        processed['geom'].append(float(matchobj.group('y')))
        processed['geom'].append(float(matchobj.group('z')))
        return ''

    nat = 0
    reconstitute = []
    processed = {}
    processed['geom'] = []
    processed['elbl'] = []

    if strict:
        for iln, line in enumerate(string.split('\n')):
            line = line.strip()
            if iln == 0:
                line = re.sub(xyz1strict, '', line)
            elif iln == 1:
                continue
            else:
                line = re.sub(atom_cartesian_strict, process_atom_cartesian, line)
            if line:
                reconstitute.append(line)
    else:
        for iln, line in enumerate(string.split('\n')):
            line = line.strip()
            if iln == 0:
                line = re.sub(xyz1, process_bohrang, line)
            elif iln == 1:
                line = re.sub(xyz2, process_system_cgmp, line)
            else:
                line = re.sub(atom_cartesian, process_atom_cartesian, line)
            if line and iln != 1:
                reconstitute.append(line)

    if 'units' not in processed:
        processed['units'] = 'Angstrom'

    #if len(processed['geom']) != nat:
    #    raise ValidationError
    #processed['fragment_files'] = []  # no EFP  # shouldn't be needed
    #processed['hint_types'] = []  # no EFP  # shouldn't be needed
    processed['geom_hints'] = []  # no EFP

    return '\n'.join(reconstitute), processed
