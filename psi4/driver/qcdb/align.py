#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import time
import itertools
import collections

import numpy as np

import qcelemental as qcel

from .molecule import Molecule
from .util import *
from .psiutil import *


class AlignmentMill(collections.namedtuple('AlignmentMill', 'shift rotation atommap mirror')):
    """Facilitates the application of the simple transformation operations
    defined by namedtuple of arrays as recipe to the data structures
    describing Cartesian molecular coordinates. Attaches functions to
    transform the geometry, element list, gradient, etc. to the
    AlignmentRecipe. When `mirror` attribute (defaults to False) active,
    then molecular system can be substantively changed by procedure.

    """

    def __new__(cls, shift, rotation, atommap, mirror=False):
        return super(AlignmentMill, cls).__new__(cls, shift, rotation, atommap, mirror)

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('AlignmentMill', width=width))
        if label:
            text.append('{:^{width}}'.format(label))
        text.append('-' * width)
        text.append('Mirror:   {}'.format(self.mirror))
        text.append('Atom Map: {}'.format(self.atommap))
        text.append('Shift:    {}'.format(self.shift))
        text.append('Rotation:')
        text.append('{}'.format(self.rotation))
        text.append('-' * width)
        return ('\n'.join(text))

    def align_coordinates(self, geom, reverse=False):
        """suitable for geometry or displaced geometry"""

        algeom = np.copy(geom)
        if reverse:
            algeom = algeom.dot(self.rotation)
            algeom = algeom + self.shift
            if self.mirror:
                algeom[:, 1] *= -1.
        else:
            if self.mirror:
                algeom[:, 1] *= -1.
            algeom = algeom - self.shift
            algeom = algeom.dot(self.rotation)
        algeom = algeom[self.atommap, :]
        # mirror-wise, rsm/msr == rms/msr

        return algeom

    def align_atoms(self, ats):
        """suitable for masses, symbols, Zs, etc."""

        return ats[self.atommap]

    def align_vector(self, vec):
        """suitable for vector attached to molecule"""

        # sensible? TODO
        #alvec = np.copy(vec)
        #if self.mirror:
        #    alvec[:, 1] *= -1
        return vec.dot(self.rotation)

    def align_gradient(self, grad):
        """suitable for vector system attached to atoms"""

        # sensible? TODO
        #algrad = np.copy(grad)
        #if self.mirror:
        #    algrad[:, 1] *= -1
        algrad = grad.dot(self.rotation)
        algrad = algrad[self.atommap]

        return algrad

    def align_hessian(self, hess):
        blocked_hess = qcdb.util.blockwise_expand(hess, (3, 3), False)
        alhess = np.zeros_like(blocked_hess)

        nat = blocked_hess.shape[0]
        for iat in range(nat):
            for jat in range(nat):
                alhess[iat, jat] = (self.rotation.T).dot(blocked_hess[iat, jat].dot(self.rotation))

        alhess = alhess[np.ix_(self.atommap, self.atommap)]

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

    def align_mini_system(self, geom, uniq, reverse=False):
        """For AlignmentRecipe `ar`, apply its translation, rotation, and atom map."""

        nugeom = self.align_coordinates(geom, reverse=reverse)
        nuuniq = self.align_atoms(uniq)

        return nugeom, nuuniq


def _nre(Z, geom):
    """Nuclear repulsion energy"""

    nre = 0.
    for at1 in range(geom.shape[0]):
        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += Z[at1] * Z[at2] / dist
    return nre


def _pseudo_nre(Zhash, geom):
    """Pseudo nuclear repulsion energy where non-physical Z contrived from `Zhash`."""

    Zidx = list(set(sorted(Zhash)))
    pZ = [Zidx.index(z) for z in Zhash]
    return _nre(pZ, geom)


#         1         2         3         4         5         6         7         8         9         10


def B787(cgeom,
         rgeom,
         cuniq,
         runiq,
         do_plot=False,
         verbose=1,
         atoms_map=False,
         run_resorting=False,
         mols_align=False,
         run_to_completion=False,
         uno_cutoff=1.e-3,
         run_mirror=False):
    """Use Kabsch algorithm to find best alignment of geometry `cgeom` onto
    `rgeom` while sampling atom mappings restricted by `runiq` and `cuniq`.

    Parameters
    ----------
    rgeom : ndarray of float
        (nat, 3) array of reference/target/unchanged geometry. Assumed [a0]
        for RMSD purposes.
    cgeom : ndarray of float
        (nat, 3) array of concern/changeable geometry. Assumed [a0] for RMSD
        purposes. Must have same nat, units, and atom content as rgeom.
    runiq : ndarray of str
        (nat,) array indicating which rows (atoms) in `rgeom` are shuffleable
        without changing the molecule. Generally hashes of element symbol and
        mass are used, but could be as simple as ['C', 'H', 'H', 'D', 'H'] for
        monodeuterated methane.
    cuniq : ndarray of str
        (nat,) array indicating which rows (atoms) in `cgeom` are shuffleable.
        See `runiq` for more details. Strings and count in `cuniq` must match
        `runiq`. That is, `sorted(cuniq) == sorted(runiq)`.
    do_plot : bool, optional
        Pops up a mpl plot showing before, after, and ref geometries.
    verbose : int, optional
        Quantity of printing. 0 to silence.
    atoms_map : bool, optional
        Whether atom1 of rgeom already corresponds to atom1 of cgeom and so on.
        If `True`, no resorting will be run, parameters `runiq` and `cuniq`
        may be passed as `None`, and much time will be saved.
    run_resorting : bool, optional
        Run the resorting machinery even if unnecessary because `atoms_map=True`.
    mols_align : bool or float, optional
        Whether ref_mol and concern_mol have identical geometries by eye
        (barring orientation or atom mapping) and expected final RMSD = 0.
        If `True`, procedure is truncated when RMSD condition met, saving time.
        If float, convcrit at which search for minimium truncates.
    run_to_completion : bool, optional
        Run reorderings to completion (past RMSD = 0) even if unnecessary because
        `mols_align=True`. Used to test worst-case timings.
    uno_cutoff : float, optional
        TODO
    run_mirror : bool, optional
        Run alternate geometries potentially allowing best match to `rgeom`
        from mirror image of `cgeom`. Only run if system confirmed to
        be nonsuperimposable upon mirror reflection.

    Returns
    -------
    float, tuple
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is a AlignmentMill namedtuple with fields
        (shift, rotation, atommap, mirror) that prescribe the transformation
        from `cgeom` and the optimally aligned geometry.

    """
    # validation
    if rgeom.shape != cgeom.shape or rgeom.shape[1] != 3:
        raise ValidationError("""natom doesn't match: {} != {}""".format(rgeom.shape, cgeom.shape))
    nat = rgeom.shape[0]
    if atoms_map and runiq is None and cuniq is None:
        runiq = np.array([''] * nat)
        cuniq = np.array([''] * nat)
    if sorted(runiq) != sorted(cuniq):
        raise ValidationError("""atom subclasses unequal:\n  {}\n  {}""".format(runiq, cuniq))

    if run_mirror:
        # use aligner to check if system and its (xz-plane) mirror image are
        #   superimposible and hence whether its worth doubling the number of Kabsch
        #   runs below to check for mirror-image matches
        mcgeom = np.copy(cgeom)
        mcgeom[:, 1] *= -1.
        exact = 1.e-6
        mrmsd, msolution = B787(
            mcgeom,
            cgeom,
            cuniq,
            cuniq,
            do_plot=False,
            verbose=0,
            atoms_map=False,
            mols_align=exact,
            run_mirror=False,
            uno_cutoff=0.1)
        superimposable = mrmsd < exact
        if verbose >= 1 and superimposable:
            print(
                'Not testing for mirror-image matches (despite `run_mirror`) since system and its mirror are superimposable'
            )

    # initialization
    best_rmsd = 100.  # [A]
    ocount = 0
    hold_solution = None
    run_resorting = run_resorting or not atoms_map
    if mols_align is True:
        a_convergence = 1.e-3
    elif mols_align is False:
        a_convergence = 0.0
    else:
        a_convergence = mols_align

    # initial presentation
    atomfmt1 = """  {:6} {:16.8f} {:16.8f} {:16.8f}    {m:16.8f} {hsh:}"""
    atomfmt2 = """  {} {:16.8f} {:16.8f} {:16.8f}"""

    if verbose >= 2:
        print('<<<  Reference:')
        for at, hsh in enumerate(runiq):
            print(atomfmt2.format(runiq[at][:6], *rgeom[at]))

        print('<<<  Concern:')
        for at, hsh in enumerate(cuniq):
            print(atomfmt2.format(cuniq[at][:6], *cgeom[at]))

    # start_rmsd is nonsense if not atoms_map
    start_rmsd = np.linalg.norm(cgeom - rgeom) * qcel.constants.bohr2angstroms / np.sqrt(nat)
    if verbose >= 1:
        print('Start RMSD = {:8.4f} [A] (naive)'.format(start_rmsd))

    def _plausible_atom_orderings_wrapper(runiq, cuniq, rgeom, cgeom, run_resorting, verbose=1, uno_cutoff=1.e-3):
        """Wrapper to _plausible_atom_orderings that bypasses it (`run_resorting=False`) when
        atoms of R & C known to be ordered. Easier to put logic here because _plausible is generator.

        """
        if run_resorting:
            return _plausible_atom_orderings(runiq, cuniq, rgeom, cgeom, verbose=verbose, uno_cutoff=uno_cutoff)
        else:
            return [np.arange(rgeom.shape[0])]

    t0 = time.time()
    tc = 0.
    for ordering in _plausible_atom_orderings_wrapper(
            runiq, cuniq, rgeom, cgeom, run_resorting, verbose=verbose, uno_cutoff=uno_cutoff):
        t1 = time.time()
        ocount += 1
        npordd = np.asarray(ordering)
        orrmsd, RR, TT = kabsch_align(rgeom, cgeom[npordd, :], weight=None)
        orrmsd = np.around(orrmsd, decimals=8)

        temp_solution = AlignmentMill(TT, RR, npordd, False)
        tgeom = temp_solution.align_coordinates(cgeom, reverse=False)
        if verbose >= 4:
            print('temp geom diff\n', tgeom - rgeom)
        temp_rmsd = np.linalg.norm(tgeom - rgeom) * qcel.constants.bohr2angstroms / np.sqrt(rgeom.shape[0])
        temp_rmsd = np.around(temp_rmsd, decimals=8)
        t2 = time.time()
        tc += t2 - t1

        if temp_rmsd < best_rmsd:
            best_rmsd = temp_rmsd
            hold_solution = temp_solution
            if verbose >= 1:
                print('<<<  trial {:8}  {} yields RMSD {}  >>>'.format(ocount, npordd, temp_rmsd))
            if not run_to_completion and best_rmsd < a_convergence:
                break
        else:
            if verbose >= 3:
                print('     trial {:8}  {} yields RMSD {}'.format(ocount, npordd, temp_rmsd))

        if run_mirror and not superimposable:
            t1 = time.time()
            ocount += 1
            icgeom = np.copy(cgeom)
            icgeom[:, 1] *= -1.
            orrmsd, RR, TT = kabsch_align(rgeom, icgeom[npordd, :], weight=None)
            orrmsd = np.around(orrmsd, decimals=8)

            temp_solution = AlignmentMill(TT, RR, npordd, True)
            tgeom = temp_solution.align_coordinates(cgeom, reverse=False)
            if verbose >= 4:
                print('temp geom diff\n', tgeom - rgeom)
            temp_rmsd = np.linalg.norm(tgeom - rgeom) * qcel.constants.bohr2angstroms / np.sqrt(rgeom.shape[0])
            temp_rmsd = np.around(temp_rmsd, decimals=8)
            t2 = time.time()
            tc += t2 - t1

            if temp_rmsd < best_rmsd:
                best_rmsd = temp_rmsd
                hold_solution = temp_solution
                if verbose >= 1:
                    print('<<<  trial {:8}m {} yields RMSD {}  >>>'.format(ocount - 1, npordd, temp_rmsd))
                if not run_to_completion and best_rmsd < a_convergence:
                    break
            else:
                if verbose >= 3:
                    print('     trial {:8}m {} yields RMSD {}'.format(ocount - 1, npordd, temp_rmsd))

    t3 = time.time()
    if verbose >= 1:
        print('Total time [s] for {:6} iterations: {:.3}'.format(ocount, t3 - t0))
        print('Hungarian time [s] for atom ordering: {:.3}'.format(t3 - t0 - tc))
        print('Kabsch time [s] for mol alignment:    {:.3}'.format(tc))

    ageom, auniq = hold_solution.align_mini_system(cgeom, cuniq, reverse=False)
    final_rmsd = np.linalg.norm(ageom - rgeom) * qcel.constants.bohr2angstroms / np.sqrt(nat)
    assert (abs(best_rmsd - final_rmsd) < 1.e-3)

    if verbose >= 1:
        print('Final RMSD = {:8.4f} [A]'.format(final_rmsd))
        print('Mirror match:', hold_solution.mirror)
        print(hold_solution)

    # final presentation & plotting
    if verbose >= 2:
        print('<<<  Aligned:')
        for at, hsh in enumerate(auniq):
            print(atomfmt2.format(auniq[at][:6], *ageom[at]))
        print('<<<  Aligned Diff:')
        for at, hsh in enumerate(auniq):
            print(atomfmt2.format(auniq[at][:6], * [ageom[at][i] - rgeom[at][i] for i in range(3)]))

    if do_plot:
        plot_coord(ref=rgeom, cand=ageom, orig=cgeom, comment='Final RMSD = {:8.4f}'.format(final_rmsd))

    # sanity checks
    compare_values(
        _pseudo_nre(cuniq, cgeom), _pseudo_nre(auniq, ageom), 4, 'D: concern_mol-->returned_mol pNRE uncorrupted', verbose=verbose-1)
    if mols_align is True:
        compare_values(
            _pseudo_nre(runiq, rgeom),
            _pseudo_nre(auniq, ageom), 4, 'D: concern_mol-->returned_mol pNRE matches ref_mol', verbose=verbose-1)
        compare_integers(True, np.allclose(rgeom, ageom, atol=4),
                         'D: concern_mol-->returned_mol geometry matches ref_mol', verbose=verbose-1)
        compare_values(0., final_rmsd, 4, 'D: null RMSD', verbose=verbose-1)

    return final_rmsd, hold_solution


def _plausible_atom_orderings(ref, current, rgeom, cgeom, algo='hunguno', verbose=1, uno_cutoff=1.e-3):
    """

    Parameters
    ----------
    ref : list
        Hashes encoding distinguishable non-coord characteristics of reference
        molecule. Namely, atomic symbol, mass, basis sets?.
    current : list
        Hashes encoding distinguishable non-coord characteristics of

    Returns
    -------
    iterator of tuples

    """
    if sorted(ref) != sorted(current):
        raise ValidationError("""ref and current can't map to each other.\n""" + 'R:  ' + str(ref) + '\nC:  ' +
                              str(current))

    where = collections.defaultdict(list)
    for iuq, uq in enumerate(ref):
        where[uq].append(iuq)

    cwhere = collections.defaultdict(list)
    for iuq, uq in enumerate(current):
        cwhere[uq].append(iuq)

    connect = collections.OrderedDict()
    for k in where:
        connect[tuple(where[k])] = tuple(cwhere[k])

    def filter_permutative(rgp, cgp):
        """Original atom ordering generator for like subset of atoms (e.g., all carbons).
        Relies on permutation. Filtering depends on similarity of structure (see `atol` parameter).
        Only suitable for total system size up to about 20 atoms.

        """
        if verbose >= 1:
            print("""Space:     {} <--> {}""".format(rgp, cgp))
        bnbn = [rrdistmat[first, second] for first, second in zip(rgp, rgp[1:])]
        for pm in itertools.permutations(cgp):
            cncn = [ccdistmat[first, second] for first, second in zip(pm, pm[1:])]
            if np.allclose(bnbn, cncn, atol=0.5):
                if verbose >= 1:
                    print('Candidate:', rgp, '<--', pm)
                yield pm

    def filter_hungarian(rgp, cgp):
        if verbose >= 1:
            print("""Space:     {} <--> {}""".format(rgp, cgp))
        submatCR = crdistmat[np.ix_(cgp, rgp)]  # this one gets manipulated by hungarian call
        submatCRcopy = np.copy(submatCR)
        lapCR = hungarian.lap(submatCR)
        ptsCR = list(zip(lapCR[1], range(len(rgp))))
        sumCR = sum(submatCRcopy[lapCR[1], range(len(rgp))])

        subans = lapCR[1]
        npcgp = np.array(cgp)
        ans = tuple(npcgp[np.array(subans)])
        if verbose >= 1:
            print('Best Candidate ({:6.3}):'.format(sumCR), rgp, '<--', ans, '     from', cgp, subans)
        yield ans

    def filter_hungarian_uno(rgp, cgp):
        """Hungarian algorithm on cost matrix based off headless (all Z same w/i space anyways) NRE.
        Having found _a_ solution and the reduced cost matrix, this still isn't likely to produce
        atom rearrangement fit for Kabsch b/c internal coordinate cost matrix doesn't nail down
        distance-equivalent atoms with different Cartesian coordinates like Cartesian-distance-matrix
        cost matrix does. So, form a bipartite graph from all essentially-zero connections between
        ref and concern and run Uno algorithm to enumerate them.

        """
        if verbose >= 1:
            print("""Space:     {} <--> {}""".format(rgp, cgp))

        # formulate cost matrix from internal (not Cartesian) layouts of R & C
        npcgp = np.array(cgp)
        submatCC = ccnremat[np.ix_(cgp, cgp)]
        submatRR = rrnremat[np.ix_(rgp, rgp)]
        sumCC = 100. * np.sum(submatCC, axis=0)  # cost mat small if not scaled, this way like Z=Neon
        sumRR = 100. * np.sum(submatRR, axis=0)
        cost = np.zeros((len(cgp), len(rgp)))
        for j in range(cost.shape[1]):
            for i in range(cost.shape[0]):
                cost[i, j] = (sumCC[i] - sumRR[j])**2
        if verbose >= 2:
            print('Cost:\n', cost)
        costcopy = np.copy(cost)  # other one gets manipulated by hungarian call

        # find _a_ best match btwn R & C atoms through Kuhn-Munkres (Hungarian) algorithm
        t00 = time.time()
        # exactly like `scipy.optimize.linear_sum_assignment(cost)` only with extra return
        #    import scipy.optimize
        #    raise ImportError("""Python module scipy >=0.17 not found. Solve by installing it: `conda install scipy` or `pip install scipy`""")
        (row_ind, col_ind), reducedcost = linear_sum_assignment(cost, return_cost=True)
        ptsCR = list(zip(row_ind, col_ind))
        ptsCR = sorted(ptsCR, key=lambda tup: tup[1])
        sumCR = costcopy[row_ind, col_ind].sum()
        t01 = time.time()
        if verbose >= 2:
            print('Reduced cost:\n', cost)
        if verbose >= 1:
            print('Hungarian time [s] for space:         {:.3}'.format(t01 - t00))

        # final _all_ best matches btwn R & C atoms through Uno algorithm, seeded from Hungarian sol'n
        edges = np.argwhere(reducedcost < uno_cutoff)
        gooduns = uno(edges, ptsCR)
        t02 = time.time()
        if verbose >= 1:
            print('Uno time [s] for space:               {:.3}'.format(t02 - t01))

        for gu in gooduns:
            gu2 = gu[:]
            gu2.sort(key=lambda x: x[1])  # resorts match into (r, c) = (info, range)
            subans = [p[0] for p in gu2]  # compacted to subans/lap format

            ans = tuple(npcgp[np.array(subans)])
            if verbose >= 3:
                print('Best Candidate ({:6.3}):'.format(sumCR), rgp, '<--', ans, '     from', cgp, subans)
            yield ans

    if algo == 'perm':
        ccdistmat = qcel.util.distance_matrix(cgeom, cgeom)
        rrdistmat = qcel.util.distance_matrix(rgeom, rgeom)
        algofn = filter_permutative

    if algo == 'hung':
        crdistmat = qcel.util.distance_matrix(cgeom, rgeom)
        algofn = filter_hungarian

    if algo == 'hunguno':
        ccdistmat = qcel.util.distance_matrix(cgeom, cgeom)
        rrdistmat = qcel.util.distance_matrix(rgeom, rgeom)
        # TODO investigate soundness
        with np.errstate(divide='ignore'):
            ccnremat = np.reciprocal(ccdistmat)
            rrnremat = np.reciprocal(rrdistmat)
        ccnremat[ccnremat == np.inf] = 0.
        rrnremat[rrnremat == np.inf] = 0.
        algofn = filter_hungarian_uno
        from .util.scipy_hungarian import linear_sum_assignment
        from .util.gph_uno_bipartite import uno

    # collect candidate atom orderings from algofn for each of the atom classes,
    #   recombine the classes with each other in every permutation (could maybe
    #   add Hungarian here, too) as generator back to permutation_kabsch
    for cpmut in itertools.product(*itertools.starmap(algofn, connect.items())):
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
    rgeom : ndarray of float
        (nat, 3) array of reference/target/unchanged geometry. Assumed [a0]
        for RMSD purposes.
    cgeom : ndarray of float
        (nat, 3) array of concern/changeable geometry. Assumed [a0] for RMSD
        purposes. Must have same Natom, units, and 1-to-1 atom ordering as rgeom.
    weight : ndarray of float
        (nat,) array of weights applied to `rgeom`. Note that definitions of
        weights (nothing to do with atom masses) are several, and I haven't
        seen one yet that can make centroid the center-of-mass and
        also make the RMSD match the usual mass-wtd-RMSD definition.
        Also, only one weight vector used rather than split btwn R & C,
        which may be invalid if not 1-to-1. Weighting is not recommended.

    Returns
    -------
    float, ndarray, ndarray
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is (3, 3) rotation matrix to optimal alignment.
        Third item is (3,) translation vector [a0] to optimal alignment.

    Sources
    -------
    Kabsch: Acta Cryst. (1978). A34, 827-828 http://journals.iucr.org/a/issues/1978/05/00/a15629/a15629.pdf
    C++ affine code: https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
    weighted RMSD: http://www.amber.utah.edu/AMBER-workshop/London-2015/tutorial1/
    protein wRMSD code: https://pharmacy.umich.edu/sites/default/files/global_wrmsd_v8.3.py.txt
    quaternion: https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures

    Author: DAS

    """
    if weight is None:
        w = np.ones((rgeom.shape[0]))
    elif isinstance(weight, list):
        w = np.asarray(weight)
    elif isinstance(weight, np.ndarray):
        w = weight
    else:
        raise ValidationError("""Unrecognized argument type {} for kwarg 'weight'.""".format(type(weight)))

    R = rgeom
    C = cgeom
    N = rgeom.shape[0]
    if np.allclose(R, C):
        # can hit a mixed non-identity translation/rotation, so head off
        return 0., np.identity(3), np.zeros(3)

    Rcentroid = R.sum(axis=0) / N
    Ccentroid = C.sum(axis=0) / N
    R = np.subtract(R, Rcentroid)
    C = np.subtract(C, Ccentroid)

    R *= np.sqrt(w[:, None])
    C *= np.sqrt(w[:, None])

    RR = kabsch_quaternion(C.T, R.T)
    TT = Ccentroid - RR.dot(Rcentroid)

    C = C.dot(RR)
    rmsd = np.linalg.norm(R - C) * qcel.constants.bohr2angstroms / np.sqrt(np.sum(w))

    return rmsd, RR, TT


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
    F = np.zeros((4, 4))
    # diagonal
    F[0, 0] = cov[0, 0] + cov[1, 1] + cov[2, 2]
    F[1, 1] = cov[0, 0] - cov[1, 1] - cov[2, 2]
    F[2, 2] = -cov[0, 0] + cov[1, 1] - cov[2, 2]
    F[3, 3] = -cov[0, 0] - cov[1, 1] + cov[2, 2]
    # Upper & lower triangle
    F[1, 0] = F[0, 1] = cov[1, 2] - cov[2, 1]
    F[2, 0] = F[0, 2] = cov[2, 0] - cov[0, 2]
    F[3, 0] = F[0, 3] = cov[0, 1] - cov[1, 0]
    F[2, 1] = F[1, 2] = cov[0, 1] + cov[1, 0]
    F[3, 1] = F[1, 3] = cov[0, 2] + cov[2, 0]
    F[3, 2] = F[2, 3] = cov[1, 2] + cov[2, 1]

    # Compute ew, ev of F
    ew, ev = np.linalg.eigh(F)

    # Construct optimal rotation matrix from leading ev
    q = ev[:, -1]
    U = np.zeros((3, 3))

    U[0, 0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    U[0, 1] = 2 * (q[1] * q[2] - q[0] * q[3])
    U[0, 2] = 2 * (q[1] * q[3] + q[0] * q[2])
    U[1, 0] = 2 * (q[1] * q[2] + q[0] * q[3])
    U[1, 1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    U[1, 2] = 2 * (q[2] * q[3] - q[0] * q[1])
    U[2, 0] = 2 * (q[1] * q[3] - q[0] * q[2])
    U[2, 1] = 2 * (q[2] * q[3] + q[0] * q[1])
    U[2, 2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2

    return U


def compute_scramble(nat, do_resort=True, do_shift=True, do_rotate=True, deflection=1.0, do_mirror=False):
    """Generate a random or directed translation, rotation, and atom shuffling.

    Parameters
    ----------
    nat : int
        Number of atoms for which to prepare an atom mapping.
    do_resort : bool or array-like, optional
        Whether to randomly shuffle atoms (`True`) or leave 1st atom 1st, etc. (`False`)
        or shuffle according to specified (nat, ) indices (e.g., [2, 1, 0])
    do_shift : bool or array-like, optional
        Whether to generate a random atom shift on interval [-3, 3) in each
        dimension (`True`) or leave at current origin (`False`) or shift along
        specified (3, ) vector (e.g., np.array([0., 1., -1.])).
    do_rotate : bool or array-like, optional
        Whether to generate a random 3D rotation according to algorithm of Arvo (`True`)
        or leave at current orientation (`False`) or rotate with specified (3, 3) matrix.
    deflection : float, optional
        If `do_rotate`, how random a rotation: 0.0 is no change, 0.1 is small
        perturbation, 1.0 is completely random.
    do_mirror : bool, optional
        Whether to set mirror reflection instruction. Changes identity of
        molecule so off by default.

    Returns
    -------
    tuple
        AlignmentMill namedtuple with fields (shift, rotation, atommap, mirror)
        as requested: identity, random, or specified.

    """
    rand_elord = np.arange(nat)
    if do_resort is True:
        np.random.shuffle(rand_elord)
    elif do_resort is False:
        pass
    else:
        rand_elord = np.array(do_resort)
        assert (rand_elord.shape == (nat, ))

    if do_shift is True:
        rand_shift = 6 * np.random.random_sample((3, )) - 3
    elif do_shift is False:
        rand_shift = np.zeros((3, ))
    else:
        rand_shift = np.array(do_shift)
        assert (rand_shift.shape == (3, ))

    if do_rotate is True:
        rand_rot3d = random_rotation_matrix(deflection=deflection)
    elif do_rotate is False:
        rand_rot3d = np.identity(3)
    else:
        rand_rot3d = np.array(do_rotate)
        assert (rand_rot3d.shape == (3, 3))

    perturbation = AlignmentMill(rand_shift, rand_rot3d, rand_elord, do_mirror)
    return perturbation
