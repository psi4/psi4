from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import collections
import numpy as np
from .molecule import Molecule
from .util import *
from .psiutil import *

try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass


AlignmentRecipe = collections.namedtuple('AlignmentRecipe', 'shift rotation atommap')

class AlignmentMill(AlignmentRecipe):

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('AlignmentMill', width=width))
        if label:
            text.append('{:^{width}}'.format(label))
        text.append('-' * width)
        text.append('Atom Map: {}'.format(self.atommap))
        text.append('Shift:    {}'.format(self.shift))
        text.append('Rotation:')
        text.append('{}'.format(self.rotation))
        text.append('-' * width)
        return('\n'.join(text))

    def align_coordinates(self, geom, reverse=False):
        """suitable for geometry or displaced geometry"""

        if reverse:
            algeom = geom.dot(self.rotation)
            algeom = algeom + self.shift
        else:
            algeom = geom - self.shift
            algeom = algeom.dot(self.rotation)
        algeom = algeom[self.atommap, :]

        return algeom

    def align_atoms(self, ats):
        """suitable for masses, symbols, Zs, etc."""

        return ats[self.atommap]

    def align_vector(self, vec):
        """suitable for vector attached to molecule"""

        return vec.dot(self.rotation)

    def align_gradient(self, grad):
        """suitable for vector system attached to atoms"""

        algrad = grad.dot(self.rotation)
        algrad = algrad[self.atommap]

        return algrad

    def align_hessian(self, hess):
        # UNTESTED

        blocked_hess = qcdb.util.blockwise_expand(hess, (3, 3), False)
        alhess = np.zeros_like(blocked_hess)

        nat = blocked_hess.shape[0]
        for iat in range(nat):
            for jat in range(nat):
                alhess[iat, jat] = (self.rotation.T).dot(blocked_hess[iat, jat].dot(self.rotation))

        alhess = alhess[:, self.atommap]
        alhess = alhess[self.atommap, :]

        alhess = qcdb.util.blockwise_contract(alhess)
        return alhess

    def align_system(self, geom, mass, elem, elez, uniq, reverse=False):
        """For AlignmentRecipe `ar`, apply its translation, rotation, and atom map."""
        
        nugeom = self.align_coordinates(geom, reverse=reverse)
        numass = self.align_atoms(mass)
        nuelem = self.align_atoms(elem)
        nuelez = self.align_atoms(elez)
        nuuniq = self.align_atoms(uniq)
        
        return nugeom, numass, nuelem, nuelez, nuuniq


def _nre(Z, geom):
    """Nuclear repulsion energy"""

    nre = 0.
    for at1 in range(geom.shape[0]):
        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += Z[at1] * Z[at2] / dist
    return nre


def _mol2pieces(mol):
    import hashlib

    mol.update_geometry()
    print(mol.nuclear_repulsion_energy())
    geom = np.asarray(mol.geometry())
    mass = np.asarray([mol.mass(at) for at in range(mol.natom())])
    elem = np.asarray([mol.symbol(at) for at in range(mol.natom())])
    elez = np.asarray([mol.Z(at) for at in range(mol.natom())])
    uniq = np.asarray([hashlib.sha1((str(elem[at]) + str(mass[at])).encode('utf-8')).hexdigest() for at in range(mol.natom())])

    return geom, mass, elem, elez, uniq


def _pieces2mol(geom, mass, elem, elez):
    """Crude, chgmult-less, fragment-less qcdb.Molecule"""

    mol = Molecule()
    mol.lock_frame = False
    mol.PYmove_to_com = False
    mol.PYfix_orientation = True

    for at in range(geom.shape[0]):
        mol.add_atom(elez[at], geom[at][0], geom[at][1], geom[at][2], elem[at], mass[at], elez[at])

    mol.set_molecular_charge(0)
    mol.set_multiplicity(1)
    mol.fragments.append([0, geom.shape[0] - 1])
    mol.fragment_types.append('Real')
    mol.fragment_charges.append(0)
    mol.fragment_multiplicities.append(1)
    mol.input_units_to_au = 1.0
    mol.set_units('Bohr')

    mol.update_geometry()
    return mol


#         1         2         3         4         5         6         7         8         9         10


def B787(ref_mol, concern_mol, atoms_map=False, mols_align=False, do_plot=False):
    """Finds shift, rotation, and atom reordering of `concern_mol` that best
    aligns with `ref_mol`.

    Parameters
    ----------
    atoms_map : bool, optional
        Whether atom1 of ref_mol corresponds to atom1 of concern_mol, etc.
        If true, specifying `True` can save much time.
    mols_align : bool, optional
        Whether ref_mol and concern_mol have identical geometries by eye
        (barring orientation or atom mapping) and expected final RMSD = 0.
    do_plot : bool, optional
        Pops up a mpl plot showing before, after, and ref geometries.

    Returns
    -------
    float, tuple, qcdb.Molecule
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is a AlignmentMill namedtuple with fields
        (shift, rotation, atommap) that prescribe the transformation from
        `concern_mol` and the optimally aligned geometry.
        Third item is a crude charge-, multiplicity-, fragment-less Molecule
        at optimally aligned (and atom-ordered) geometry.

    """
    import collections
    from psi4.driver.qcdb.physconst import psi_bohr2angstroms

    atomfmt = """  {:6} {:16.8f} {:16.8f} {:16.8f}    {m:16.8f} {hsh:}"""

    rgeom, rmass, relem, relez, runiq = _mol2pieces(ref_mol)
    print('<<<  Reference:', _nre(relez, rgeom))
    for at, sym in enumerate(relem):
        print(atomfmt.format(sym, *rgeom[at], m=rmass[at], hsh=runiq[at]))

    cgeom, cmass, celem, celez, cuniq = _mol2pieces(concern_mol)
    print('<<<  Concern:', _nre(celez, cgeom))
    for at, sym in enumerate(celem):
        print(atomfmt.format(sym, *cgeom[at], m=cmass[at], hsh=cuniq[at]))

    if rgeom.shape != cgeom.shape:
        raise ValidationError("""natom doesn't match""")
    nat = rgeom.shape[0]

    rmsd = np.linalg.norm(cgeom - rgeom) * psi_bohr2angstroms / np.sqrt(nat)
    print('Start RMSD = {:8.4f}'.format(rmsd))

    if not atoms_map and nat > 20:
        raise ValidationError('Stop right there. Permutation algorithm is not ready for systems this large.')

    ret_rmsd, solution = permutative_kabsch(rgeom, cgeom, runiq, cuniq, atoms_map=atoms_map)
    print(ret_rmsd)
    print(solution)

    ageom, amass, aelem, aelez, auniq = solution.align_system(cgeom, cmass, celem, celez, cuniq, reverse=False)
    amol = _pieces2mol(ageom, amass, aelem, aelez)

    print('<<<  Aligned:', _nre(aelez, ageom))
    for at, sym in enumerate(aelem):
        print(atomfmt.format(sym, *ageom[at], m=amass[at], hsh=auniq[at]))

    rmsd = np.linalg.norm(ageom - rgeom) * psi_bohr2angstroms / np.sqrt(nat)
    if do_plot:
        plot_coord(ref=rgeom, cand=ageom, orig=cgeom, comment='Final RMSD = {:8.4f}'.format(rmsd))

    compare_values(_nre(celez, cgeom), _nre(aelez, ageom), 4, 'D: concern_mol-->returned_mol NRE uncorrupted')
    compare_values(concern_mol.nuclear_repulsion_energy(), amol.nuclear_repulsion_energy(), 4, 'Q: concern_mol-->returned_mol NRE uncorrupted')
    if mols_align:
        compare_values(_nre(relez, rgeom), _nre(aelez, ageom), 4, 'D: concern_mol-->returned_mol NRE matches ref_mol')
        compare_values(ref_mol.nuclear_repulsion_energy(), amol.nuclear_repulsion_energy(), 4, 'Q: concern_mol-->returned_mol NRE matches ref_mol')
        compare_integers(True, np.allclose(rgeom, ageom, atol=4), 'D: concern_mol-->returned_mol geometry matches ref_mol')
        compare_integers(True, np.allclose(np.asarray(ref_mol.geometry()), np.asarray(amol.geometry()), atol=4), 'Q: concern_mol-->returned_mol geometry matches ref_mol')

    return rmsd, solution, amol


def _plausible_atom_rearrangements(ref, current, rgeom, cgeom):
    """

    Parameters
    ----------
    ref idhash : list
        Hashes encoding distinguishable non-coord characteristics of molecule.
        Namely, atomic symbol, mass, basis sets?.

    Returns
    -------
    iterator of tuples

    """
    import itertools
    import collections

    try:
        assert(sorted(ref) == sorted(current))
    except AssertionError:
        raise qcdb.ValidationError("""ref and current can't map to each other.\n""" + 'R:  ' + str(ref) + '\nC:  ' + str(current))

    where = collections.defaultdict(list)
    for iuq, uq in enumerate(ref):
        where[uq].append(iuq)

    cwhere = collections.defaultdict(list)
    for iuq, uq in enumerate(current):
        cwhere[uq].append(iuq)

    connect = collections.OrderedDict()
    for k in where:
        connect[tuple(where[k])] = tuple(cwhere[k])

    for k, v in connect.items():
        print('ASDF', k, v)
#        print("""{:30} <--> {:30}""".format(k, v))

    def distance_matrix(a, b):
        assert(a.shape[1] == b.shape[1])
        distm = np.zeros([a.shape[0], b.shape[0]])
        for i in range(a.shape[0]):
            for j in range(b.shape[0]):
                distm[i, j] = np.linalg.norm(a[i] - b[j])
        return distm

    def filter_more_geom(rgp, cgp):
        for pm in itertools.permutations(cgp):
            bnbn = [rrdistmat[first, second] for first, second in zip(rgp, rgp[1:])]
            cncn = [ccdistmat[first, second] for first, second in zip(pm, pm[1:])]
            if np.allclose(bnbn, cncn, atol=0.5):
                print('Candidate:', rgp, '<--', pm)
                yield pm

    # ccdistmat = scipy.spatial.distance.cdist(cgeom, cgeom, 'euclidean')
    # rrdistmat = scipy.spatial.distance.cdist(rgeom, rgeom, 'euclidean')
    ccdistmat = distance_matrix(cgeom, cgeom)
    rrdistmat = distance_matrix(rgeom, rgeom)

    for cpmut in itertools.product(*itertools.starmap(filter_more_geom, connect.items())):
        atpat = [None] * len(ref)
        for igp, group in enumerate(cpmut):
            for iidx, idx in enumerate(list(connect.keys())[igp]):
                atpat[idx] = group[iidx]
        yield atpat


def kabsch_align(rgeom, cgeom, weight=None):
    """Finds optimal translation and rotation to align `cgeom` onto `rgeom` via
    Kabsch algorithm by minimizing the norm of the residual, || R - U * C ||.

    Parameters
    ----------
    rgeom : np.array(float)
        Natom x 3 array of reference/target/unchanged geometry. Assumed [a0]
        for RMSD purposes.
    cgeom : np.array(float)
        Natom x 3 array of concern/changeable geometry. Assumed [a0] for RMSD
        purposes. Must have same Natom, units, and 1-to-1 atom ordering as rgeom.
    weight : np.array(float)
        Natom array of weights applied to `rgeom`. Note that definitions of
        weights (nothing to do with atom masses) are several, and I haven't
        seen one yet that can make centroid the center-of-mass and
        also make the RMSD match the usual mass-wtd-RMSD definition.
        Also, only one weight vector used rather than split btwn R & C,
        which may be invalid if not 1-to-1. Weighting is not recommended.

    Returns
    -------
    float, np.array, np.array
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is 3 x 3 rotation matrix to optimal alignment.
        Third item is 1 x 3 translation vector [a0] to optimal alignment.

    Sources
    -------
    Kabsch: Acta Cryst. (1978). A34, 827-828 http://journals.iucr.org/a/issues/1978/05/00/a15629/a15629.pdf
    C++ affine code: https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
    weighted RMSD: http://www.amber.utah.edu/AMBER-workshop/London-2015/tutorial1/
    protein wRMSD code: https://pharmacy.umich.edu/sites/default/files/global_wrmsd_v8.3.py.txt
    quaternion: https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures

    Author: DAS

    """
    from psi4.driver.qcdb.physconst import psi_bohr2angstroms

    if weight is None:
        w = np.ones((rgeom.shape[0]))
    elif isinstance(weight, list):
        w = np.asarray(weight)
    elif isinstance(weight, np.ndarray):
        w = weight
    else:
        raise qcdb.ValidationError(
            """Unrecognized argument type {} for kwarg 'weight'.""".format(type(weight)))

    R = rgeom
    C = cgeom
    N = rgeom.shape[0]

    Rcentroid = R.sum(axis=0) / N
    Ccentroid = C.sum(axis=0) / N
    R = R - Rcentroid
    C = C - Ccentroid

    R = R * np.sqrt(w[:, None])
    C = C * np.sqrt(w[:, None])

    RR = kabsch_quaternion(C.T, R.T)
    TT = Ccentroid - RR.dot(Rcentroid)

    C = C.dot(RR)
    rmsd = np.linalg.norm(R - C) * psi_bohr2angstroms / np.sqrt(np.sum(w))

    return rmsd, RR, TT


def permutative_kabsch(rgeom, cgeom, runiq, cuniq, atoms_map=False):
    """Use Kabsch algorithm to find best alignment of geometry `cgeom` onto
    `rgeom` while sampling atom mappings restricted by `runiq` and `cuniq`.

    Parameters
    ----------
    rgeom : np.array(float)
        Natom x 3 array of reference/target/unchanged geometry. Assumed [a0]
        for RMSD purposes.
    cgeom : np.array(float)
        Natom x 3 array of concern/changeable geometry. Assumed [a0] for RMSD
        purposes. Must have same Natom, units, and atom content as rgeom.
    runiq : np.array(str)
        Natom array indicating which rows (atoms) in `rgeom` are shuffleable
        without changing the molecule. Generally hashes of element symbol and
        mass are used, but could be as simple as ['C', 'H', 'H', 'D', 'H'] for
        monodeuterated methane.
    cuniq : np.array(str)
        Natom array indicating which rows (atoms) in `cgeom` are shuffleable.
        See `runiq` for more details. Strings and count in `cuniq` must match
        `runiq`. That is, `sorted(cuniq) == sorted(runiq)`.
    atoms_map : bool, optional
        Whether atom1 of rgeom already corresponds to atom1 of cgeom and so on.
        If `True`, no permutation will be run, parameters `runiq` and `cuniq`
        may be passed as `None`, and much time will be saved.

    Returns
    -------
    float, tuple
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is a AlignmentMill namedtuple with fields
        (shift, rotation, atommap) that prescribe the transformation from
        `cgeom` and the optimally aligned geometry.

    """
    import time
    from psi4.driver.qcdb.physconst import psi_bohr2angstroms

    t0 = time.time()
    best_rmsd = 100.  # [A]
    ocount = 0
    hold_solution = None

    def atomrr(runiq, cuniq, rgeom, cgeom, trivial):
        if trivial:
            return [np.arange(rgeom.shape[0])]
        else:
            return _plausible_atom_rearrangements(runiq, cuniq, rgeom, cgeom)

    for ordd in atomrr(runiq, cuniq, rgeom, cgeom, trivial=atoms_map):
        ocount += 1
        npordd = np.asarray(ordd)
        orrmsd, RR, TT = kabsch_align(rgeom, cgeom[npordd, :], weight=None)
        orrmsd = np.around(orrmsd, decimals=8)

        temp_solution = AlignmentMill(TT, RR, npordd)
        #tgeom = align_coordinates(temp_solution, cgeom, reverse=False)
        tgeom = temp_solution.align_coordinates(cgeom, reverse=False)

        temp_rmsd = np.linalg.norm(tgeom - rgeom) * psi_bohr2angstroms / np.sqrt(rgeom.shape[0])
        temp_rmsd = np.around(temp_rmsd, decimals=8)
    #         print('---  trial {} {} yields RMSD {} K {} L <? {}'.format(ocount, npordd, orrmsd, temp_rmsd, best_rmsd))

        if temp_rmsd < best_rmsd:
            best_rmsd = temp_rmsd
            hold_solution = temp_solution
            print('<<<  trial {} {} yields RMSD {}'.format(ocount, npordd, temp_rmsd))
    #         if best_rmsd < 0.2:
    #             break

    t2 = time.time()
    print('Kabsch + Permute time', ocount, 'is', t2-t0)

    return best_rmsd, hold_solution


def kabsch_quaternion(P, Q):
    """Computes the optimal rotation matrix U which mapping a set of points P
    onto the set of points Q according to the minimization of || Q - R * P ||,
    using the unit quaternion formulation of the Kabsch algorithm.

    Arguments:
    <np.ndarray> P := MxN array. M=dimension of space, N=number of points.
    <np.ndarray> Q := MxN array. M=dimension of space, N=number of points.

    Returns:
    <np.ndarray> U := Optimal MxM rotation matrix mapping P onto Q.

    Author: DAS

    """
    # Form covariance matrix
    cov = Q.dot(P.T)

    # Form the quaternion transformation matrix F
    F = np.zeros((4,4))
    # diagonal
    F[0,0] = cov[0,0] + cov[1,1] + cov[2,2]
    F[1,1] = cov[0,0] - cov[1,1] - cov[2,2]
    F[2,2] = -cov[0,0] + cov[1,1] - cov[2,2]
    F[3,3] = -cov[0,0] - cov[1,1] + cov[2,2]
    # Upper & lower triangle
    F[1,0] = F[0,1] = cov[1,2] - cov[2,1]
    F[2,0] = F[0,2] = cov[2,0] - cov[0,2]
    F[3,0] = F[0,3] = cov[0,1] - cov[1,0]
    F[2,1] = F[1,2] = cov[0,1] + cov[1,0]
    F[3,1] = F[1,3] = cov[0,2] + cov[2,0]
    F[3,2] = F[2,3] = cov[1,2] + cov[2,1]

    # Compute ew, ev of F
    ew, ev = np.linalg.eigh(F)

    # Construct optimal rotation matrix from leading ev
    q = ev[:,-1]
    U = np.zeros((3,3))

    U[0,0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    U[0,1] = 2*(q[1]*q[2] - q[0]*q[3])
    U[0,2] = 2*(q[1]*q[3] + q[0]*q[2])
    U[1,0] = 2*(q[1]*q[2] + q[0]*q[3])
    U[1,1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    U[1,2] = 2*(q[2]*q[3] - q[0]*q[1])
    U[2,0] = 2*(q[1]*q[3] - q[0]*q[2])
    U[2,1] = 2*(q[2]*q[3] + q[0]*q[1])
    U[2,2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2

    return U
