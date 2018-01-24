import re

import numpy as np

from ..exceptions import *


def from_arrays(geom,
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
                fragment_types=None,
                fragment_charges=None,
                fragment_multiplicities=None,
                molecular_charge=None,
                molecular_multiplicity=None,
                nonphysical=False,
                mtol=1.e-3,
                verbose=1):
    """Compose a Molecule dict from unvalidated arrays and variables, returning dict.

    minimum is geom and one of elem, elez, elbl


    Parameters
    ----------
    See fields of return  molrec below. Required parameters are `geom` and one of `elem`, `elez`, `elbl`
    geom : array-like
        (nat, 3) or (3 * nat, ) ndarray or list o'lists of Cartesian coordinates.
    fragment_separators : array-like of int, optional
        (nfr - 1, ) list of atom indices at which to split `geom` into fragments.
    fragment_types : array-like of {'Real', 'Ghost', 'Absent'}, optional
        (nfr, ) list of fragment types. If given, size must be consistent with `fragment_separators`.
    elbl : ndarray of str
        (nat, ) Label extending `elem` symbol, possibly conveying isotope, mass, tagging information.

    nonphysical : bool, optional


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
        If `units='Angstrom'`, overrides consumer's value for [Å]-->[a0] conversion.
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
    fragment_types : list of {'Real', 'Ghost', 'Absent'}
        (nfr, ) list of fragment types.
    fragment_charges : list of float
        (nfr, ) list of charge allocated to each fragment.
    fragment_multiplicities : list of int
        (nfr, ) list of multiplicity allocated to each fragment.
    molecular_charge : float
        total charge on system.
    molecular_multiplicity : int
        total multiplicity on system.

    """
    from .chgmult import validate_and_fill_chgmult

    # TODO not handled: elbl, fragment_*, validation of chgmult overall vs frag

    molinit = {}
    molinit.update(validate_and_fill_universals(name=name,
                                                units=units,
                                                input_units_to_au=input_units_to_au,
                                                fix_com=fix_com,
                                                fix_orientation=fix_orientation,
                                                fix_symmetry=fix_symmetry))

    molinit.update(validate_and_fill_geometry(  geom=geom))
    nat = molinit['geom'].shape[0] // 3

    molinit.update(validate_and_fill_nuclei(    nat,
                                                elea=elea,
                                                elez=elez,
                                                elem=elem,
                                                mass=mass,
                                                real=real,
                                                elbl=elbl,
                                                nonphysical=nonphysical,
                                                mtol=mtol,
                                                verbose=verbose))

    molinit.update(validate_and_fill_fragments( nat,
                                                fragment_separators=fragment_separators,
                                                fragment_types=fragment_types,
                                                fragment_charges=fragment_charges,
                                                fragment_multiplicities=fragment_multiplicities))

    Z_available = molinit['elez'] * molinit['real'] * 1.
    molinit.update(validate_and_fill_chgmult(   zeff=Z_available,
                                                fragment_separators=molinit['fragment_separators'],
                                                molecular_charge=molecular_charge,
                                                fragment_charges=molinit['fragment_charges'],
                                                molecular_multiplicity=molecular_multiplicity,
                                                fragment_multiplicities=molinit['fragment_multiplicities'],
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


def validate_and_fill_geometry(geom=None):
    """Check `geom` for overlapping atoms. Return flattened"""

    from ..util import distance_matrix

    npgeom = np.array(geom, dtype=np.float).reshape((-1, 3))
    dm = distance_matrix(npgeom, npgeom)

    iu = np.triu_indices(dm.shape[0])
    dm[iu] = 10.
    tooclose = np.where(dm < 0.1)

    if tooclose[0].shape[0]:
        raise ValidationError("""Following atoms are too close: {}""".format([(i, j, trial[i, j]) for i, j in zip(*tooclose)]))

    return {'geom': npgeom.reshape((-1))}


def validate_and_fill_nuclei(nat,
                             elea=None,
                             elez=None,
                             elem=None,
                             mass=None,
                             real=None,
                             elbl=None,
                             # processing details
                             nonphysical=False,
                             mtol=1.e-3,
                             verbose=1):
    """Check the nuclear identity arrays for consistency and fill in knowable values."""

    from .nucleus import reconcile_nucleus

    if elea is None:
        elea = np.full((nat,), None)
    else:
        elea = np.array(elea)

    if elez is None:
        elez = np.full((nat,), None)
    else:
        elez = np.array(elez)

    if elem is None:
        elem = np.full((nat,), None)
    else:
        elem = np.array(elem)

    if mass is None:
        mass = np.full((nat,), None)
    else:
        mass = np.array(mass)

    if real is None:
        real = np.full((nat,), None)
    else:
        real = np.array(real)

    if elbl is None:
        elbl = np.full((nat,), None)
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
                """Fragment quantities given without separation info: sep ({}), types ({}), chg ({}), and mult ({})""".format(
                fragment_separators, fragment_types, fragment_charges, fragment_multiplicities))
    else:
        trial_geom = np.zeros((nat, 3))
        try:
            split_geom = np.split(trial_geom, fragment_separators, axis=0)
        except TypeError:
            raise ValidationError(
                """fragment_separators ({}) unable to perform trial np.split on geometry.""".format(fragment_separators))
        if any(len(f) == 0 for f in split_geom):
            raise ValidationError(
                """fragment_separators ({}) yields zero-length fragment(s) after trial np.split on geometry.""".format(split_geom))
        if sum(len(f) for f in split_geom) != nat:
            raise ValidationError(
                """fragment_separators ({}) yields overlapping fragment(s) after trial np.split on geometry, possibly unsorted.""".format(split_geom))
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
        elif all(f is None or (isinstance(f, (int, np.int64)) and f >= 1) for f in fragment_multiplicities):
            frm = fragment_multiplicities
        else:
            raise ValidationError("""fragment_multiplicities not among None or positive integer: {}""".format(fragment_multiplicities))

    if not (len(frt) == len(frc) == len(frm) == len(frs) + 1):
        raise ValidationError("""Dimension mismatch among fragment quantities: sep + 1 ({}), types ({}), chg ({}), and mult({})""".format(
            frs.shape[0] + 1, len(frt), len(frc), len(frm)))

    return {'fragment_separators': list(frs),
            'fragment_types': frt,
            'fragment_charges': frc,
            'fragment_multiplicities': frm}
